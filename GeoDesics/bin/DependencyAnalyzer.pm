package DependencyAnalyzer;

use strict;
use warnings FATAL => 'all';

use Cwd;
use File::Path;
use IPC::Open2;
use IPC::Open3;

use Exporter qw(import);
our @EXPORT_OK = qw(demangle get_baseclass library_sort parse_makefile);

# On Macs with Xcode CLT version >= 7.3 the nm output is incorrectly formatted,
# so we prefer to use the 'nm-classic' executable if it is present
# NOTE: nm-classic has no version argument, so use `type` to test its existence
my $MacNM = '/Library/Developer/CommandLineTools/usr/bin/nm-classic';
my $NM = (system("type $MacNM >/dev/null 2>&1")==0) ? $MacNM : 'nm';

# Dependenxy database version number.  Increment whenever the format
# of any of the files changes to force a full regeneration.
my $DBVERSION = 0;

sub new {
  my($class, $CODE_HOME, $datadir) = @_;
  my $self = {
              CODE_HOME => $CODE_HOME,
              datadir => $datadir || 'MakefileRules/.dependencies',

              cxx => '',
              debug_symbols => undef,

              spec_libraries => undef,

              symbols => undef,
              deps => undef,
              factory => undef,

              symbol_files => undef,
              library_deps => undef,

              symbol_cache => {},
             };
  if(!$self->{CODE_HOME}) {
    chomp($self->{CODE_HOME} = `git rev-parse --show-toplevel`);
  }
  $self->{CANONICAL_CODE_HOME} = Cwd::abs_path($self->{CODE_HOME});
  mkpath("$self->{CANONICAL_CODE_HOME}/$self->{datadir}");
  bless $self, $class;
}

# Accessors/mutators

sub cxx {
  my $self = shift;
  if(@_) {
    $self->{cxx} = join ' ', @_;
  }
  return $self->{cxx};
}

sub debug_symbol {
  my $self = shift;
  @{$self->{debug_symbols}}{@_} = ();
}

# map file => [ symbols ]
sub symbols {
  my $self = shift;
  if(!defined $self->{symbols}) {
    if(defined $self->{symbol_files}) {
      foreach my $symbol (sort keys %{$self->{symbol_files}}) {
        push @{$self->{symbols}->{$self->{symbol_files}->{$symbol}}}, $symbol;
      }
    } else {
      $self->{symbols} = read_map($self->filename('symbols'));
    }
  }
  return wantarray ? %{$self->{symbols}} : $self->{symbols};
}

# map file => [ dependency files ]
sub deps {
  my $self = shift;
  $self->{deps} = read_map($self->filename('deps'))
    unless defined $self->{deps};
  return wantarray ? %{$self->{deps}} : $self->{deps};
}

# map library => { file => [ symbols ] }
sub factory {
  my $self = shift;
  if(!defined $self->{factory}) {
    my $tmp = read_map($self->filename('factory'));
    while(my($lib, $obj_and_symbols_array) = each %$tmp) {
      foreach my $obj_and_symbols (@$obj_and_symbols_array) {
        my @symbols = split /\|/, $obj_and_symbols;
        my $object = shift @symbols;
        $self->{factory}->{$lib}->{$object} = \@symbols;
      }
    }
  }
  return wantarray ? %{$self->{factory}} : $self->{factory};
}

# map symbol => file
sub symbol_files {
  my $self = shift;
  local $_;
  if(!defined $self->{symbol_files}) {
    my %symbols = $self->symbols;
    $self->{symbol_files} = {};
    while(my($file, $libs) = each %symbols) {
      $self->{symbol_files}->{$_} = $file foreach @$libs;
    }
  }
  return wantarray ? %{$self->{symbol_files}} : $self->{symbol_files};
}

# dependency graph between libraries (for dynamic linking)
sub library_deps {
  my $self = shift;
  if(!defined $self->{library_deps}) {
    # collapse the dependency graph to the library level
    local $_;
    my %library_deps;
    my %object_deps = $self->deps;
    while(my($obj,$deps) = each %object_deps) {
      (my $lib = $obj) =~ s/\[.*//;
      foreach my $dep (@$deps) {
        (my $deplib = $dep) =~ s/\[.*//;
        $library_deps{$lib}{$deplib} = 1;
      }
    }
    $_ = [ sort keys %$_ ] foreach values %library_deps;
    $self->{library_deps} = \%library_deps;
  }
  return wantarray ? %{$self->{library_deps}} : $self->{library_deps};
}


# write database to disk
sub write {
  my $self = shift;
  $self->write_symbols
    if defined $self->{symbols} or defined $self->{symbol_files};
  $self->write_deps    if defined $self->{deps};
  $self->write_factory if defined $self->{factory};
}

sub write_symbols {
  my $self = shift;
  write_map(scalar $self->symbols, $self->filename('symbols'));
}

sub write_deps {
  my $self = shift;
  write_map(scalar $self->deps, $self->filename('deps'));
}

sub write_factory {
  my $self = shift;
  my %factory = $self->factory;
  my %tmp;
  while(my($lib, $hash) = each %factory) {
    while(my($obj, $symbols) = each %$hash) {
      push @{$tmp{$lib}}, join('|', $obj, @$symbols);
    }
  }
  write_map(\%tmp, $self->filename('factory'));
}

sub dependency_info_exists {
  my $self = shift;
  return(
         -e $self->filename('symbols') and
         -e $self->filename('deps') and
         -e $self->filename('factory')
        );
}

# More complicated stuff

# Convert paths relative to CODE_HOME to paths relative to here
sub prefix_paths {
  my $self = shift;
  local $_;
  if($self->{CODE_HOME} eq '.') {
    return @_;
  } else {
    return map { "$self->{CODE_HOME}/$_" } @_;
  }
}

# Main worker routine in --generate mode.  Scans all the libraries.
sub generate {
  my $self = shift;

  # If the database version has changed since the database was
  # written, nuke the whole thing and regenerate it.
  my $writtenversion = 0;
  my $versionfile = $self->filename('version');
  if(-e $versionfile) {
    open my $vfh, "<$versionfile" or die "Can't read $versionfile: $!\n";
    $writtenversion = <$vfh>;
    close $vfh;
  }
  if($writtenversion != $DBVERSION) {
    my $dbglob = $self->filename('*');
    foreach my $dbfile (<$dbglob>) {
      # Don't delete the version file so that the delete will be
      # finished next time if we are interrupted.
      next if $dbfile eq $versionfile;
      rmtree($dbfile);
    }
  }
  {
    open my $vfh, ">$versionfile" or die "Can't write $versionfile: $!\n";
    print $vfh $DBVERSION, "\n";
    close $vfh;
  }

  # The force file will force regeneration of the dependency graph if the
  # generation is interrupted, since there's no way to know whether it is in
  # sync with the symbol cache in that case.  The cache updates should be safe,
  # since the madelib files are not removed until after the cache files are
  # written.
  my $forcefile = $self->filename('force');
  my $force = (-e $forcefile) || !$self->dependency_info_exists;
  open my $touch, ">$forcefile" or die "Can't create $forcefile: $!\n";
  close $touch;

  if($self->update_symbol_caches or $force) {
    $self->process_libraries;
    $self->build_depgraph;
  }

  unlink $forcefile;
}

# updates symbol list and factory list
sub process_libraries {
  my $self = shift;

  my %symbol_defs;
  my %createhelpers;
  foreach my $lib ($self->spec_libraries) {
    my %libsymbols = %{$self->get_library_symbols($lib)->{defined}};
    while(my($file, $symbols) = each %libsymbols) {
      foreach my $symbol (@$symbols) {
        if(defined $self->{debug_symbols} and
           exists $self->{debug_symbols}->{$symbol}) {
          print $symbol, " defined in ", $file, "\n";
        }

        if(exists $symbol_defs{$symbol}) {
          $symbol_defs{$symbol} = 'MULTIPLY_DEFINED';
        } else {
          $symbol_defs{$symbol} = $file;
        }

        if($symbol =~ /registered_/) {
          $createhelpers{$file} = [] unless exists $createhelpers{$file};
        } elsif($symbol =~ /Factory.*CreateHelper/) {
          push @{$createhelpers{$file}}, $symbol;
        }
      }
    }
  }

  # We used to remove symbols that matched those in the system libraries, to
  # avoid linking in the wrong symbols. According to Will, this was only
  # necessary due to libMemTest (which no longer exists) declaring operator
  # new and being linked instead of the system operator new.
  $self->{symbol_files} = \%symbol_defs;
  undef $self->{symbols};


  my %factory_map;
  my %factory_directories;
  while(my($file,$symbols) = each %createhelpers) {
    $file =~ m#(.*)/lib([^/[]+)\.a\[#
      or die "Can't parse library from $file\n";

    # check for uniqueness
    if(exists $factory_directories{$2}) {
      if($factory_directories{$2} ne $1) {
        die "Multiple directories named $2: $1 and $factory_directories{$2}\n";
      }
    } else {
      $factory_directories{$2} = $1;
    }

    $factory_map{$2}{$file} = $symbols;
  }
  $self->{factory} = \%factory_map;
}

# updates deps
sub build_depgraph {
  my $self = shift;

  my %symbol_defs = $self->symbol_files;

  my %object_depgraph;
  foreach my $lib ($self->spec_libraries) {
    my %lib_symbols = %{$self->get_library_symbols($lib)->{undefined}};
    while(my($file, $symbols) = each %lib_symbols) {
      my %needs;
      foreach my $symbol (@$symbols) {
        if(exists $symbol_defs{$symbol}) {
          $self->die_multiply_defined($symbol,$file)
            if $symbol_defs{$symbol} eq 'MULTIPLY_DEFINED';
          $needs{$symbol_defs{$symbol}} = 1;
        }
      }
      $object_depgraph{$file} = [ sort keys %needs ] if keys %needs;
    }
  }
  $self->{deps} = \%object_depgraph;
}

# returns a dpendency graph only containing nodes used for linking $object with
# @factory_libs
sub restricted_deps {
  my $self = shift;
  my($object, @factory_libs) = @_;
  local $_;

  my %symbol_files = $self->symbol_files;
  my %object_depgraph = $self->deps;

  my %needed_objects;
  my %initial_objects;
  foreach my $symbol ($self->get_object_symbols($object, 'undefined')) {
    next unless exists $symbol_files{$symbol};
    my $def = $symbol_files{$symbol};
    $self->die_multiply_defined($symbol) if $def eq 'MULTIPLY_DEFINED';
    add_to_hash_with_dependencies(\%needed_objects, $def, \%object_depgraph);
    $initial_objects{$def} = 1;
  }
  if(@factory_libs) {
    my %spec_libs;
    my %factory_map = $self->factory;
    foreach my $lib (@factory_libs) {
      if(exists $factory_map{$lib}) {
        foreach my $file (keys %{$factory_map{$lib}}) {
          add_to_hash_with_dependencies(\%needed_objects, $file,
                                        \%object_depgraph);
          $initial_objects{$file} = 1;
        }
      } else {
        if(!keys %spec_libs) {
          foreach my $obj (values %symbol_files) {
            if($obj =~ /lib(.*)\.a\[.*/) {
              $spec_libs{$1} = 1;
            }
          }
        }
        if(exists $spec_libs{$lib}) {
          warn "No Factory objects found in $lib\n";
        } else {
          warn "Library $lib does not exist.\n";
        }
      }
    }
  }

  foreach (keys %object_depgraph) {
    delete $object_depgraph{$_} unless exists $needed_objects{$_};
  }
  foreach (keys %needed_objects) {
    $object_depgraph{$_} = [] unless exists $object_depgraph{$_};
  }

  return
    wantarray ? (\%object_depgraph, sort keys %initial_objects)
      : \%object_depgraph;
}


# Read all the external symbols from an object, and return a list of symbols.
# The second argument should be be 'defined' or 'undefined', and will cause
# only symbols matching that description to be returned.
sub get_object_symbols {
  my $self = shift;
  my($file, $which) = @_;
  if($which eq 'undefined') {
    my %symbols = %{$self->get_library_symbols($file)->{undefined}};
    die "$file is not an object\n" if keys %symbols > 1;
    return () if keys %symbols == 0;
    return @{(values %symbols)[0]};
  } elsif($which eq 'defined') {
    my %symbols = %{$self->get_library_symbols($file)->{defined}};
    my %weaksymbols = %{$self->get_library_symbols($file)->{weak_undefined}};
    die "$file is not an object\n"
      if keys %symbols > 1 or keys %weaksymbols > 1;
    my @symbols;
    push @symbols, @{(values %symbols)[0]} unless keys %symbols == 0;
    push @symbols, @{(values %weaksymbols)[0]} unless keys %weaksymbols == 0;
    return @symbols;
  }
  die "Symbol type should be 'defined' or 'undefined', not '$which'\n";
}

# Print out a report on a needed multiply defined symbol and abort.
sub die_multiply_defined {
  my $self = shift;
  my($symbol, $neededby) = @_;
  warn "Multiply defined symbol '$symbol' [".demangle($symbol)."] needed"
    . (defined $neededby?" by $neededby":'') . "\n";
  $self->debug_symbol($symbol);
  $self->process_libraries;
  die "Aborting\n";
}


# internal methods

sub filename {
  my $self = shift;
  my $name = shift;
  return "$self->{CODE_HOME}/$self->{datadir}/$name";
}

# Find all lib*.a files in the tree
sub spec_libraries {
  my $self = shift;
  if(!defined $self->{spec_libraries}) {
    local $_;

    # Silence make jobserver warnings
    # In perl 5.12 the next two lines can be 'delete local $ENV{MAKEFLAGS};'
    local %ENV = %ENV;
    delete $ENV{MAKEFLAGS};

    open LIBS, '-|', 'make', '-sC', $self->{CANONICAL_CODE_HOME},
      'allprintlibrarynames'
        or die "Can't run make allprintlibrarynames: $!\n";
    my @libs = <LIBS>;
    close LIBS or die "Can't find SpEC libraries\n";
    chomp @libs;
    $self->{spec_libraries} = [ grep { -r } @libs ];
  }
  return wantarray ? @{$self->{spec_libraries}} : $self->{spec_libraries};
}

# Update symbol caches as needed
# Returns number of libraries with ABI changes.  In particular, returns false
# if no updates were necessary.
sub update_symbol_caches {
  my $self = shift;
  my $updated = 0;
  my $base = $self->filename('madelib');
  while(my $madelib = <$base.*>) {
    (my $lib = $madelib) =~ s/^\Q$base.\E//;
    $lib =~ s#\^#/#g;
    $updated++ if $self->generate_cache($lib);
    unlink $madelib;
  }
  return $updated;
}

# Read all the external symbols from an object or library, and return a hash
# of (defined,undefined,weak_undefined) => { object => [ symbols ] }.
sub get_library_symbols {
  my $self = shift;
  my($file) = @_;
  if(!exists $self->{symbol_cache}->{$file}) {
    if($file =~ /\.(?:a|so)$/) {
      $self->{symbol_cache}->{$file} = $self->read_symbol_cache($file);
      $self->generate_cache($file)
        unless defined $self->{symbol_cache}->{$file};
    } else {
      $self->{symbol_cache}->{$file} = $self->read_symbols($file);
    }
  }

  return $self->{symbol_cache}->{$file};
}

# Generate symbol caches for a given library.
# Returns true if the cache was generated, false if doing so was unneccesary.
sub generate_cache {
  my $self = shift;
  my($file) = @_;

  my $cachedata = $self->read_symbol_cache($file);
  my $newdata = $self->read_symbols($file);

  $self->{symbol_cache}->{$file} = $newdata;

  if(symbols_differ($cachedata, $newdata)) {
    $self->write_symbol_cache($file, $newdata);
    return 1;
  } else {
    return 0;
  }
}

# Read symbol caches for the given file and return them.  Return undef if they
# do not exist.
sub read_symbol_cache {
  my $self = shift;
  my($file) = @_;

  my %symbols;
  my $slashlessfile = Cwd::abs_path($file);
  return undef unless defined $slashlessfile;
  $slashlessfile =~ s#\Q$self->{CANONICAL_CODE_HOME}\E/##;
  $slashlessfile =~ s#/#^#g;
  foreach my $type (qw(defined undefined weak_undefined)) {
    my $cachefile = $self->filename("$type.$slashlessfile");
    return undef unless -r $cachefile;
    $symbols{$type} = read_map($cachefile);
  }
  return \%symbols;
}

# Write out cached data for a file.
sub write_symbol_cache {
  my $self = shift;
  my($file,$symbols) = @_;

  my $slashlessfile = Cwd::abs_path($file);
  return undef unless defined $slashlessfile;
  $slashlessfile =~ s#\Q$self->{CANONICAL_CODE_HOME}\E/##;
  $slashlessfile =~ s#/#^#g;
  foreach my $type (qw(defined undefined weak_undefined)) {
    write_map($symbols->{$type},$self->filename("$type.$slashlessfile"));
  }
}

# Read symbol data from a file and return it.
sub read_symbols {
  my $self = shift;
  my($file) = @_;

  local $_;
  my @args = qw(-P -A -g -p);
  push @args, '--dynamic' if $file =~ /\.so(?:\.\d+)*$/;

  open my $olderr, ">&STDERR" or die "Can't dup STDERR: $!\n";
  open STDERR, '>/dev/null' or die "Can't redirect STDERR: $!\n";
  open my $nm, '-|', $NM, @args, $file or die "Can't run $NM: $!\n";
  open STDERR, '>&', $olderr or die "Can't restore STDERR\n";
  my %symbols = ( defined => {}, undefined => {}, weak_undefined => {} );
  while(<$nm>) {
    chomp;
    /^(.+): (\S+) (\S)/ or die "Don't understand nm output '$_'\n";
    my($symbol,$object,$type) = ($2,$1,$3);

    # Ignore symbols containing strange characters.
    next if $symbol =~ /\W/;

    $object =~ s#^\Q$self->{CANONICAL_CODE_HOME}\E/##;
    my $defined = $type eq 'U' ? 'undefined'
      : $type =~ /[wv]/ ? 'weak_undefined' : 'defined';
    push @{$symbols{$defined}{$object}}, $symbol;
  }
  close $nm; # Ok if this fails

  # sort so we can compare to chaced values easily
  foreach my $type (keys %symbols) {
    foreach my $object (keys %{$symbols{$type}}) {
      $symbols{$type}{$object} = [ sort @{$symbols{$type}{$object}} ];
    }
  }

  return \%symbols;
}

##### Utility functions #####

# public utility functions

# library_sort -- sorts the libraries given the dependencies among the objects
# and the initially pulled objects.
{
  my %obj_deps;
  my %obj_lib;
  my %obj_addable;
  my %lib_objs;
  my %lib_unaddable_count;
  my @sorted;

  sub library_sort {
    my($deps, @initial) = @_;
    local $_;

    %obj_deps = %$deps;
    undef %obj_lib;
    undef %obj_addable;
    undef %lib_objs;
    undef %lib_unaddable_count;
    undef @sorted;

    foreach my $obj (keys %obj_deps) {
      $obj =~ m#([^[]+)# or die "Unexpected object format: $obj\n";
      my $lib = $1;
      $lib_objs{$lib}{$obj} = 1;
      $lib_unaddable_count{$lib}++;
      $obj_lib{$obj} = $lib;
    }

    eval {
      library_sort_obj_set_addable($_) foreach @initial;
      # If we get here, we have a true circular dependency.  Add things
      # randomly until it goes away.
    circular: {
        foreach my $lib (sort keys %lib_objs) {
          foreach my $obj (sort keys %{$lib_objs{$lib}}) {
            if($obj_addable{$obj}) {
              warn "Circular dependency: adding $lib\n"
                if $ENV{DEPENDENCY_ANALYZER_VERBOSE};
              library_sort_lib_add($lib);
              redo circular;
            }
          }
        }
      }
    };
    die $@ if $@ and $@ ne "library_sort_done\n";
    die "Ordering failed.  This can't happen.\n" if keys %lib_objs != 0;
    return @sorted;
  }

  sub library_sort_obj_set_addable {
    my($obj) = @_;
    local $_;
    return if exists $obj_addable{$obj};
    $obj_addable{$obj} = 1;
    my $lib = $obj_lib{$obj};
    if(--$lib_unaddable_count{$lib} == 0) {
      library_sort_lib_add($lib);
    } else {
      foreach (@{$obj_deps{$obj}}) {
        if($obj_lib{$_} eq $lib) {
          library_sort_obj_set_addable($_);
        }
      }
    }
  }

  sub library_sort_lib_add {
    my($lib) = @_;
    local $_;
    push @sorted, $lib;
    foreach my $obj (sort keys %{$lib_objs{$lib}}) {
      if($obj_addable{$obj}) {
        delete $lib_objs{$lib}{$obj};
        library_sort_obj_set_addable($_) foreach @{$obj_deps{$obj}};
      }
    }
    delete $lib_objs{$lib} if keys %{$lib_objs{$lib}} == 0;
    die "library_sort_done\n" if keys %lib_objs == 0;
  }
}

# Demangle c++ symbol names by passing them through c++filt
{
  my($cxxfilt_out,$cxxfilt_in);
  sub demangle {
    if(!defined $cxxfilt_out) {
      my $demangle_style_opt = $^O eq 'darwin' ? '-_' : '-n';
      # Try gc++filt first, because it is probably more up-to-date if
      # it exists
      foreach my $filt (qw(gc++filt c++filt)) {
        ### This solution works in perl 5.16, but in perl 5.8.8 open2
        ### is extremely buggy and if the command to execute does not
        ### exist it both returns success and allows its child process
        ### to escape the function call, so you end up with two copies
        ### of the script running, one of which thinks the other one
        ### is a functioning c++filt.  Although I haven't tested, the
        ### version in 5.10 and 5.12 looks to fix that bug but
        ### introduce a new one where the child prints warnings.  So
        ### instead we try to figure out what will succeed before
        ### making the open2 call.
        # eval {
        #   open2($cxxfilt_out, $cxxfilt_in, $filt, $demangle_style_opt);
        #   last;
        # };
        # undef $cxxfilt_out;
        # undef $cxxfilt_in;
        if(system("$filt --version >/dev/null 2>&1") == 0) {
          open2($cxxfilt_out, $cxxfilt_in, $filt, $demangle_style_opt);
          last;
        }
      }
      die "Can't find a working c++filt\n" unless defined $cxxfilt_out;
    }
    print $cxxfilt_in $_[0], "\n";
    my $ans = <$cxxfilt_out>;
    chomp $ans;

    # Hack around gcc/binutils version mismatches.  Old versions
    # mangle argument packs using an "II", while newer ones use an
    # "IJ".  If we failed to demangle and our string contains "IJ" try
    # replacing it with "II".
    if($ans eq $_[0] and $_[0] =~ /IJ/) {
      (my $hacked_symbol = $_[0]) =~ s/IJ/II/g;
      my $hacked_ans = demangle($hacked_symbol);
      return $hacked_ans if $hacked_ans ne $hacked_symbol;
      # Still no good.  Just fall through and return the original
      # mangled name.
    }

    # Hack around ticket #951 for variadic factory functions with gcc/intel
    my $regex = qr/RKSs[IJ](.*)EE12CreateHelper/;
    if ($ans eq $_[0] and $_[0] =~ /$regex/) {
      (my $hacked_symbol = $_[0]) =~ s/$regex/RKSs$1E12CreateHelper/;
      my $hacked_ans = demangle($hacked_symbol);
      return $hacked_ans if $hacked_ans ne $hacked_symbol;
      # Still no good.  Just fall through and return the original mangled name.
    }

    # Hack around gcc 4.9 / binutils 2.20 mismatches.
    $regex = qr/J(EE12CreateHelper)/;
    if ($ans eq $_[0] and $_[0] =~ /$regex/) {
      (my $hacked_symbol = $_[0]) =~ s/$regex/I$1/;
      my $hacked_ans = demangle($hacked_symbol);
      return $hacked_ans if $hacked_ans ne $hacked_symbol;
      # Still no good.  Just fall through and return the original mangled name.
    }

    return $ans;
  }
}

# Given a mangled name of a Factory CreateHelper, CreateDerivedClass, or
# similar, return the base class name.
sub get_baseclass {
  my($name) = @_;
  my $demangled = demangle($name);
  # Conveniently, the return type of these things is always Base*
  $demangled =~ /^(.+)\* Factory::/
    or die "$demangled doesn't look like a Factory creation helper\n";
  return $1;
}


# internal utility functions

# Read $file into $hash as a map from (first element on line) to (list of
# remaining elements).
sub read_map {
  my($file) = @_;
  local $_;
  my %map;
  open my $fh, $file or die "Can't open $file: $!\n";
  while(<$fh>) {
    chomp;
    my($head, @rest) = split;
    $map{$head} = \@rest;
  }
  close $fh;
  return \%map;
}

# inverse of read_map
sub write_map {
  my($hash, $filename) = @_;
  open my $fh, ">$filename" or die "Can't open $filename: $!\n";
  while(my($key,$val) = each %$hash) {
    print $fh join(' ', $key, @$val), "\n";
  }
  close $fh;
}

# Arranges for $file and all nodes reachable from $file in the directed graph
# $depgraph to be set to 1 in $hash.
sub add_to_hash_with_dependencies {
  my($hash, $file, $depgraph) = @_;
  return if exists $hash->{$file};
  $hash->{$file} = 1;
  return unless exists $depgraph->{$file};
  foreach my $dep (@{$depgraph->{$file}}) {
    add_to_hash_with_dependencies($hash, $dep, $depgraph);
  }
}

# Returns true if the symbol lists differ (or if both are undef).
sub symbols_differ {
  my($list1,$list2) = @_;

  return 1 unless defined $list1 and defined $list2;
  foreach my $type (qw(defined undefined weak_undefined)) {
    my $objects1 = $list1->{$type};
    my $objects2 = $list2->{$type};

    return 1 if keys %$objects1 != keys %$objects2;
    foreach my $object (keys %$objects1) {
      return 1 unless exists $objects2->{$object};
      my $symbols1 = $objects1->{$object};
      my $symbols2 = $objects2->{$object};

      return 1 unless @$symbols1 == @$symbols2;
      for(my $i=0;$i<$#$symbols1;$i++) {
        return 1 unless $symbols1->[$i] eq $symbols2->[$i];
      }
    }
  }
  return 0;
}

# Parse file lists out of a Makefile.  Returns undef on failure.
sub parse_makefile {
  my($makefile) = @_;
  open MAKEFILE, $makefile or return undef;
  my $state = '';
  my %lists;
  while(<MAKEFILE>) {
    if(/\$\(NULL\)/) {
      $state = '';
    } elsif($state ne '') {
      push @{$lists{$state}}, /\s*(.*?)\s*\\/;
    } elsif(/^(\S+)\s*=\s*\\/) {
      $state = $1;
    }
  }
  close MAKEFILE;
  return \%lists;
}

1;
