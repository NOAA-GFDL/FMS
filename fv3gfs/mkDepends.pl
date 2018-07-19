#!/usr/bin/env perl

# Modifications to Brian Eaton's original to relax the restrictions on
# source file name matching module name and only one module per source
# file.  See the new "-m" and "-d" options for details.
#
# One important limitation remains.  If your module is named "procedure",
# this script will quietly ignore it.
#
#   Tom Henderson
#   Global Systems Division, NOAA/OAR
#   Mar 2011
#
# Brian Eaton's original comments follow:
#
# Generate dependencies in a form suitable for inclusion into a Makefile.
# The source filenames are provided in a file, one per line.  Directories
# to be searched for the source files and for their dependencies are provided
# in another file, one per line.  Output is written to STDOUT.
#
# For CPP type dependencies (lines beginning with #include) the dependency
# search is recursive.  Only dependencies that are found in the specified
# directories are included.  So, for example, the standard include file
# stdio.h would not be included as a dependency unless /usr/include were
# one of the specified directories to be searched.
#
# For Fortran module USE dependencies (lines beginning with a case
# insensitive "USE", possibly preceded by whitespace) the Fortran compiler
# must be able to access the .mod file associated with the .o file that
# contains the module.  In order to correctly generate these dependencies
# two restrictions must be observed.
# 1) All modules must be contained in files that have the same base name as
#    the module, in a case insensitive sense.  This restriction implies that
#    there can only be one module per file.
# 2) All modules that are to be contained in the dependency list must be
#    contained in one of the source files in the list provided on the command
#    line.
# The reason for the second restriction is that since the makefile doesn't
# contain rules to build .mod files the dependency takes the form of the .o
# file that contains the module.  If a module is being used for which the
# source code is not available (e.g., a module from a library), then adding
# a .o dependency for that module is a mistake because make will attempt to
# build that .o file, and will fail if the source code is not available.
#
# Author: B. Eaton
#         Climate Modelling Section, NCAR
#         Feb 2001

use Getopt::Std;
use File::Basename;

# Check for usage request.
@ARGV >= 2                          or usage();

# Process command line.
my %opt = ();
getopts( "t:wmd:", \%opt )        or usage();
my $filepath_arg = shift()        or usage();
my $srcfile_arg = shift()         or usage();
@ARGV == 0                        or usage();  # Check that all args were processed.

my $obj_dir;
if ( defined $opt{'t'} ) { $obj_dir = $opt{'t'}; }

my $additional_obj = "";
if ( defined $opt{'d'} ) { $additional_obj = $opt{'d'}; }

open(FILEPATH, $filepath_arg) or die "Can't open $filepath_arg: $!\n";
open(SRCFILES, $srcfile_arg) or die "Can't open $srcfile_arg: $!\n";

# Make list of paths to use when looking for files.
# Prepend "." so search starts in current directory.  This default is for
# consistency with the way GNU Make searches for dependencies.
my @file_paths = <FILEPATH>;
close(FILEPATH);
chomp @file_paths;
unshift(@file_paths,'.');
foreach $dir (@file_paths) {  # (could check that directories exist here)
    $dir =~ s!/?\s*$!!;  # remove / and any whitespace at end of directory name
    ($dir) = glob $dir;  # Expand tildes in path names.
}

# Make list of files containing source code.
my @src = <SRCFILES>;
close(SRCFILES);
chomp @src;

my %module_files = ();

#TODO:  DRY this out
if ( defined $opt{'m'} ) {
    # Attempt to parse each file for /^\s*module/ and extract module names
    # for each file.
    my ($f, $name, $path, $suffix, $mod);
    my @suffixes = ('\.[fF]90', '\.[fF]' );
    foreach $f (@src) {
        ($name, $path, $suffix) = fileparse($f, @suffixes);
        open(FH, $f)  or die "Can't open $f: $!\n";
        while ( <FH> ) {
            # Search for module definitions.
            if ( /^\s*MODULE\s+(\w+)/i ) {
                ($mod = $1) =~ tr/a-z/A-Z/;
                # skip "module procedure foo" statements
                if ( $mod ne "PROCEDURE" ) {
                    if ( defined $module_files{$mod} ) {
                        die "Duplicate definitions of module $mod in $module_files{$mod} and $name: $!\n";
                    }
                $module_files{$mod} = $path.$name;
                }
            }
        }
        close( FH );
    }
} else {
    # For each file that may contain a Fortran module (*.[fF]90 *.[fF]) convert the
    # file's basename to uppercase and use it as a hash key whose value is the file's
    # basename.  This allows fast identification of the files that contain modules.
    # The only restriction is that the file's basename and the module name must match
    # in a case insensitive way.
    my ($f, $name, $path, $suffix, $mod);
    my @suffixes = ('\.[fF]90', '\.[fF]' );
    foreach $f (@src) {
        ($name, $path, $suffix) = fileparse($f, @suffixes);
        ($mod = $name) =~ tr/a-z/A-Z/;
        $module_files{$mod} = $path.$name;
    }
}

#print STDERR "\%module_files\n";
#while ( ($k,$v) = each %module_files ) {
#    print STDERR "$k => $v\n";
#}

# Find module and include dependencies of the source files.
my ($file_path, $rmods, $rincs);
my %file_modules = ();
my %file_includes = ();
my @check_includes = ();
foreach $f ( @src ) {

    # Find the file in the seach path (@file_paths).
    unless ($file_path = find_file($f)) {
	if (defined $opt{'w'}) {print STDERR "$f not found\n";}
	next;
    }

    # Find the module and include dependencies.
    ($rmods, $rincs) = find_dependencies( $file_path );

    # Remove redundancies (a file can contain multiple procedures that have
    # the same dependencies).
    $file_modules{$f} = rm_duplicates($rmods);
    $file_includes{$f} = rm_duplicates($rincs);

    # Make a list of all include files.
    push @check_includes, @{$file_includes{$f}};
}

#print STDERR "\%file_modules\n";
#while ( ($k,$v) = each %file_modules ) {
#    print STDERR "$k => @$v\n";
#}
#print STDERR "\%file_includes\n";
#while ( ($k,$v) = each %file_includes ) {
#    print STDERR "$k => @$v\n";
#}
#print STDERR "\@check_includes\n";
#print STDERR "@check_includes\n";

# Find include file dependencies.
my %include_depends = ();
while (@check_includes) {
    $f = shift @check_includes;
    if (defined($include_depends{$f})) { next; }

    # Mark files not in path so they can be removed from the dependency list.
    unless ($file_path = find_file($f)) {
	$include_depends{$f} = -1;
	next;
    }

    # Find include file dependencies.
    ($rmods, $include_depends{$f}) = find_dependencies($file_path);

    # Add included include files to the back of the check_includes list so
    # that their dependencies can be found.
    push @check_includes, @{$include_depends{$f}};

    # Add included modules to the include_depends list.
    if ( @$rmods ) { push @{$include_depends{$f}}, @$rmods;  }
}

#print STDERR "\%include_depends\n";
#while ( ($k,$v) = each %include_depends ) {
#    print STDERR (ref $v ? "$k => @$v\n" : "$k => $v\n");
#}

# Remove include file dependencies that are not in the Filepath.
my $i, $ii;
foreach $f (keys %include_depends) {

    unless (ref $include_depends{$f}) { next; }
    $rincs = $include_depends{$f};
    unless (@$rincs) { next; }
    $ii = 0;
    $num_incs = @$rincs;
    for ($i = 0; $i < $num_incs; ++$i) {
    	if ($include_depends{$$rincs[$ii]} == -1) {
	    splice @$rincs, $ii, 1;
	    next;
	}
    ++$ii;
    }
}

# Substitute the include file dependencies into the %file_includes lists.
foreach $f (keys %file_includes) {
    my @expand_incs = ();

    # Initialize the expanded %file_includes list.
    my $i;
    unless (@{$file_includes{$f}}) { next; }
    foreach $i (@{$file_includes{$f}}) {
	push @expand_incs, $i  unless ($include_depends{$i} == -1);
    }
    unless (@expand_incs) {
	$file_includes{$f} = [];
	next;
    }

    # Expand
    for ($i = 0; $i <= $#expand_incs; ++$i) {
	push @expand_incs, @{ $include_depends{$expand_incs[$i]} };
    }

    $file_includes{$f} = rm_duplicates(\@expand_incs);
}

#print STDERR "expanded \%file_includes\n";
#while ( ($k,$v) = each %file_includes ) {
#    print STDERR "$k => @$v\n";
#}

# Print dependencies to STDOUT.
foreach $f (sort keys %file_modules) {
    $f =~ /(.+)\./;
    $target = "$1.o";
    if ( defined $opt{'t'} ) { $target = "$opt{'t'}/$1.o"; }
    push(@{$file_modules{$f}},$additional_obj);
    print "$target : $f @{noncircular(@{$file_modules{$f}},$target)} @{noncircular(@{$file_includes{$f}},$target)}\n";
}

#--------------------------------------------------------------------------------------

sub noncircular
{
  # Return an array identical to that represented by the first argument, except
  # for the absence of the element specified by the second argument.
  my @a=();
  my $x=pop(@_);
  foreach (@_) { unless ($_ eq $x) { push(@a,$_) } };
  return \@a;
}

sub find_dependencies {

    # Find dependencies of input file.
    # Use'd Fortran 90 modules are returned in \@mods.
    # Files that are "#include"d by the cpp preprocessor are returned in \@incs.

    my( $file ) = @_;
    my( @mods, @incs );

    open(FH, $file)  or die "Can't open $file: $!\n";

    while ( <FH> ) {
	# Search for "#include" and strip filename when found.
	if ( /^#include\s+[<"](.*)[>"]/ ) {
	     push @incs, $1;
	 }
	# Search for module dependencies.
	elsif ( /^\s*USE\s+(\w+)/i ) {
	    # Return dependency in the form of a .o version of the file that contains
	    # the module.
	    ($module = $1) =~ tr/a-z/A-Z/;
	    if ( defined $module_files{$module} ) {
	        if ( defined $obj_dir ) {
		    push @mods, "$obj_dir/$module_files{$module}.o";
	        } else {
	            push @mods, "$module_files{$module}.o";
	        }
	    }
	}
    }
    close( FH );
    return (\@mods, \@incs);
}

#--------------------------------------------------------------------------------------

sub find_file {

# Search for the specified file in the list of directories in the global
# array @file_paths.  Return the first occurance found, or the null string if
# the file is not found.

    my($file) = @_;
    my($dir, $fname);

    foreach $dir (@file_paths) {
	$fname = "$dir/$file";
	if ( -f  $fname ) { return $fname; }
    }
    return '';  # file not found
}

#--------------------------------------------------------------------------------------

sub rm_duplicates {

# Return a list with duplicates removed.

    my ($in) = @_;       # input arrary reference
    my @out = ();
    my $i;
    my %h = ();
    foreach $i (@$in) {
	$h{$i} = '';
    }
    @out = keys %h;
    return \@out;
}

#--------------------------------------------------------------------------------------

sub usage {
    ($ProgName = $0) =~ s!.*/!!;            # name of program
    die <<EOF
SYNOPSIS
     $ProgName [-t dir] [-w] [-m] [-d objfile] Filepath Srcfiles
OPTIONS
     -t dir
          Target directory.  If this option is set the .o files that are
          targets in the dependency rules have the form dir/file.o.
     -w   Print warnings to STDERR about files or dependencies not found.
     -m   Search source files for module definitions instead of assuming
          that module names match file base names.  With this option, more
          than one module is allowed to reside in a single source file.
     -d objfile
          Additional object file to be added to every dependence.
ARGUMENTS
     Filepath is the name of a file containing the directories (one per
     line) to be searched for dependencies.  Srcfiles is the name of a
     file containing the names of files (one per line) for which
     dependencies will be generated.
EOF
}
