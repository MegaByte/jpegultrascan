#!/usr/bin/perl
=head1 NAME
jpegultrascan
=head1 DESCRIPTION
JPEG lossless recompressor that tries all scan possibilities to minimize size
=head1 VERSION
1.3.3 2021-03-18
=head1 LICENSE
Copyright 2015-2021 Aaron Kaluszka

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

    http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.
=cut

use strict;
use warnings;
use open ':std', IO => ':raw';
use File::Spec;
use File::Temp 'tempfile';
use Getopt::Long qw(:config bundling);

my $maxbits = 3;
my $threads = 1;
my $jpegtran = my $jpegtran2 = 'jpegtran';
my $home = glob '~';
GetOptions 'q' => \my $quiet, 'v' => \my $verbose, 'i:i' => \my $incompat, 'j' => \my $app0, 's' => \my $strip, 'a' => \my $arith,
           'p=s' => \$jpegtran, 'o=s' => \$jpegtran2, 'b=i' => \$maxbits, 't:i' => \$threads, 'h|?' => \my $help;
my @strip = ('-copy', defined $strip ? $strip eq 'some' ? 'comments' : 'none' : 'all');
$arith = defined $arith ? '-arithmetic' : '';
$verbose = defined $verbose;
$app0 = defined $app0;
$incompat = defined $incompat ? $incompat == 0 ? 1 : $incompat : 0;
$threads = 8 if $threads < 1;
$maxbits = 3 if $maxbits > 3;
$maxbits = 0 if $maxbits < 0;
$jpegtran =~s /^~/$home/;
$jpegtran2 =~s /^~/$home/;

open STDOUT, '>', File::Spec->devnull if defined $quiet;
open my $STDOUT, '>&STDOUT';

@ARGV == 2 && !defined $help or die 'usage: jpegultrascan.pl [switches] inputfile outputfile
JPEG lossless recompressor that tries all scan possibilities to minimize size
Switches:
  -s      Strip all extra markers
  -s some Strip only non-comment markers
  -i      Allow multiple planes per DC scan, which may be incompatible with some
          software
  -i -1   Disallow single scan DC components, which may be incompatible with
          some software
  -j      Strip APP0 segment for 18-byte savings (generates non-compliant JPEG
          that may be incompatible with some software)
  -a      Use arithmetic coding (unsupported by most software; jpegtran support
          required)
  -b 0-3  Maximum number of bit splits to test (default: 3)
  -t [N]  Number of simultaneous processes (default: 8 if specified,
          1 otherwise)
  -q      Suppress all output
  -v      Verbose output
  -p path Path to jpegtran (default: jpegtran)
  -o path Path to other jpegtran (default: jpegtran)
';

my ($fin, $fout) = @ARGV;
$fin =~s /^~/$home/;
$fout =~s /^~/$home/;
my $insize = -s $fin;
my (undef, $ftmp) = tempfile;
my (undef, $jtmp) = tempfile;
my $maxperfile = 60;
my $minchunk = $threads;
my (%sizes, $planes, %optimal, %optimal2, @scans, %dcsize, %dcscans, @acsize, @acscans, $acallsize, $acallscans, $width, $width2);
my ($bestsize, $best) = (~0, '');
my ($tran, $tran2);
undef $/;
$| = 1;

sub read_file($) {
  my ($file) = @_;
  open my $FILE, $file or die "Couldn't read file $file";
  my $data = <$FILE>;
  close $FILE;
  $data;
}

sub write_file($@) {
  my $file = shift;
  open my $FILE, '>', $file or die "Couldn't write file $file";
  print $FILE join '', @_;
  close $FILE;
}

sub tran(@) {
  `"$jpegtran" -optimize @_`;
}

sub wintran(@) {
  my (undef, $tmp) = tempfile SUFFIX => $$;
  system qq{"$jpegtran" -optimize @_ "$tmp"};
  my $data = read_file $tmp;
  unlink $tmp;
  $data;
}

sub scansizes($) {
  my ($file) = @_;
  my @sizes;
  my $scanstart = 0;
  while ($file =~ /\xFF+(.)/gs) {
    my $type = ord $1;
    next unless $type;
    my $len = ($type == 217 ? 0 : unpack 'n', substr $file, $+[0], 2) + 2;
    if ($scanstart) {
      $sizes[$#sizes] += $-[0] - $scanstart;
      last if $type == 217;
      push @sizes, $len;
      $scanstart = 0;
    } elsif (!@sizes && $type >= 192 && $type <= 207 && $type != 196 && $type != 200 && $type != 204) {
      push @sizes, 0;
    } elsif (@sizes) {
      $sizes[$#sizes] += $len;
    }
    $scanstart = $-[0] + $len if $type == 218;
  }
  @sizes;
}

sub transcan($$) {
  my ($start, $end) = @_;
  my (undef, $tmp) = tempfile SUFFIX => $$;
  my @sizes;
  for my $i ($start .. $end) {
    write_file $tmp, @{$scans[$i]};
    my $data = &$tran('-scans', "\"$tmp\"", $arith, "\"$jtmp\"");
    if (length $data == 0) {
      print STDERR "\nNo data returned for:\n", join '', @{$scans[$i]};
      next;
    }
    my @sizesi = scansizes $data;
    if ($#sizesi != $#{$scans[$i]}) {
      print STDERR "\nIncorrect encoding for:\n", join '', @{$scans[$i]};
      next;
    }
    if ($verbose) {
      for (0 .. $#sizesi) {
        my $scan = @{$scans[$i]}[$_];
        if (!defined $sizes{$scan}) {
          $sizes{$scan} = $sizesi[$_];
          chop $scan;
          printf $STDOUT "%-${width}s %${width2}d\n", $scan, $sizesi[$_];
        }
      }
    } else {
      print $STDOUT '.';
    }
    push @sizes, $i, @sizesi;
  }
  unlink $tmp;
  @sizes;
}

sub pipefork() {
  pipe my $parent, my $child or die "Couldn't pipe\n";
  my $pid = fork;
  die "fork failed: $!" unless defined $pid;
  if ($pid) {
    close $child;
  } else {
    close $parent;
    close STDOUT;
    open STDOUT, '>&=', fileno $child or die "Couldn't open STDOUT\n";
  }
  ($parent, $pid);
}

sub doscans() {
  if ($threads == 1 || $#scans < $minchunk) {
    my @sizes = transcan 0, $#scans;
    unless ($verbose) {
      for (my $i = 0; $i < $#sizes; ++$i) {
        $sizes{$_} = $sizes[++$i] for @{$scans[$sizes[$i]]};
      }
    }
  } else {
    my $chunk = int @scans / $threads;
    $chunk = $minchunk if $chunk < $minchunk;
    my $end = scalar @scans;
    my @fh;
    for (my $start = 0; $start < $end; $start += $chunk) {
	  my ($fh, $pid) = pipefork;
      unless ($pid) {
        print join ',', transcan $start, ($start + $chunk < $end - $minchunk ? $start + $chunk : $end) - 1;
        exit;
      }
      push @fh, $fh;
    }
    for (@fh) {
      my @sizes = split ',', <$_>;
      for (my $i = 0; $i < $#sizes; ++$i) {
        $sizes{$_} = $sizes[++$i] for @{$scans[$sizes[$i]]};
      }
    }
  }
  @scans = ();
}

sub queuescans(@) {
  push @scans, \@_;
}

sub queuescansmulti(@) {
  my (@dc, @ac);
  for (my $i = 0; $i < @_; ++$i) {
    if (@dc + @ac > $maxperfile) {
      do {
        pop @ac;
        --$i;
      } while $ac[$#ac] !~ /0;\n$/;
      queuescans @dc, @ac;
      @ac = ();
    }
    if ($_[$i] =~ /: 0 0/) {
      push @dc, splice @_, $i, 1;
      --$i;
    } else {
      push @ac, $_[$i];
    }
  }
  queuescans @dc, @ac if @ac > 0;
}

sub generatebitscans($$$$) {
  my ($plane, $a, $b, $numbits) = @_;
  my @scans = "$plane: $a $b 0 $numbits;\n";
  for (my $i = $numbits; $i > 0; --$i) {
    push @scans, "$plane: $a $b $i " . ($i - 1) . ";\n";
  }
  @scans;
}

sub generatedcscans($) {
  my ($plane) = @_;
  my $allplanes = join ' ', 0 .. $planes;
  my $dummy;
  if ($plane ne $allplanes) {
    $dummy = $allplanes;
    $dummy =~s /$_ ?// for split ' ', $plane;
    $dummy =~s /^(.*?) ?$/$1: 0 0 0 9;\n/;
  }
  for (0 .. $maxbits) {
    my @scans = generatebitscans $plane, 0, 0, $_;
    push @scans, $dummy if defined $dummy;
    queuescans @scans;
  }
}

sub generateacscans($) {
  my ($plane) = @_;
  my $dummy = join(' ', 0 .. $planes) . ": 0 0 0 9;\n";
  for (my $numcoeffs = 0; $numcoeffs < 63; ++$numcoeffs) {
    for (my $offset = 0; $offset < $numcoeffs + 1 && $offset < 63 - $numcoeffs; ++$offset) {
      for (my $numbits = 0; $numbits <= $maxbits; ++$numbits) {
        my @scans = $dummy;
        push @scans, generatebitscans $plane, 1, $offset, $numbits if $offset > 0;
        my $i = $offset + 1;
        for (; $i + $numcoeffs < 63; $i += $numcoeffs + 1) {
          push @scans, generatebitscans $plane, $i, $i + $numcoeffs, $numbits;
        }
        push @scans, generatebitscans $plane, $i, 63, $numbits if $i <= 63;
        @scans <= $maxperfile ? queuescans @scans : queuescansmulti @scans;
      }
    }
  }
}

sub choosebest($$$$) {
  my ($leftsize, $left, $rightsize, $right) = @_;
  return if !$leftsize && !$rightsize;
  return ($leftsize, $left) if !$rightsize || $leftsize && $leftsize < $rightsize;
  return ($rightsize, $right) if !$leftsize || $rightsize && $leftsize > $rightsize || scalar(() = $left =~ /\n/g) > scalar(() = $right =~ /\n/g);
  ($leftsize, $left);
}

sub trybitsplits($$$$$);
sub trybitsplits($$$$$) {
  my ($plane, $a, $b, $c, $d) = @_;
  my $key = "$plane $a $b $c $d";
  unless (defined $optimal2{$key}) {
    my $minscans = "$plane: $a $b $c $d;\n";
    my $minsize = $sizes{$minscans};
    for ($a .. $b - 1) {
      my ($asize, $ascans) = trybitsplits $plane, $a, $_, $c, $d;
      my ($bsize, $bscans) = trybitsplits $plane, $_ + 1, $b, $c, $d;
      ($minsize, $minscans) = choosebest $asize + $bsize, $ascans . $bscans, $minsize, $minscans;
    }
    $optimal2{$key} = [$minsize, $minscans];
  }
  @{$optimal2{$key}};
}

sub trysplits($$$);
sub trysplits($$$) {
  my ($plane, $a, $b) = @_;
  my $key = "$plane $a $b";
  unless (defined $optimal{$key}) {
    my ($minsize, $minscans) = (~0, '');
    for (0 .. $maxbits) {
      my ($allsize, $allscans) = trybitsplits $plane, $a, $b, 0, $_;
      for (my $i = $_; $i > 0; --$i) {
        my ($size, $scans) = trybitsplits $plane, $a, $b, $i, $i - 1;
        $allsize += $size;
        $allscans .= $scans;
      }
      ($minsize, $minscans) = choosebest $allsize, $allscans, $minsize, $minscans;
    }
    for ($a .. $b - 1) {
      my ($asize, $ascans) = trysplits $plane, $a, $_;
      my ($bsize, $bscans) = trysplits $plane, $_ + 1, $b;
      ($minsize, $minscans) = choosebest $asize + $bsize, $ascans . $bscans, $minsize, $minscans;
    }
    $optimal{$key} = [$minsize, $minscans];
  }
  @{$optimal{$key}};
}

sub app0remove($) {
  my ($file) = @_;
  while ("\xFF\xE0" eq substr $file, 2, 2) {
    my $size = 2 + unpack 'n', substr $file, 4, 2;
    substr $file, 2, $size, '';
    print $verbose ? "\nRemoved $size APP0 bytes.\n" : '.';
  }
  $file;
}

sub partitionscans($@);
sub partitionscans($@) {
  my $offset = shift;
  my ($dcallsize, $dcallscans) = (0, '');
  for (@_) {
    if (!defined $dcsize{$_}) {
      generatedcscans $_;
      doscans;
      ($dcsize{$_}, $dcscans{$_}) = trysplits $_, 0, 0;
    }
    $dcallsize += $dcsize{$_};
    $dcallscans .= $dcscans{$_};
  }
  ($bestsize, $best) = choosebest $dcallsize + $acallsize, $dcallscans . $acallscans, $bestsize, $best;
  my @all = map "$_;\n", @_;
  for (@all) {
    if (!defined $sizes{$_}) {
      queuescans @all;
      doscans;
      last;
    }
  }
  my $allsize = 0;
  for (@all) {
    $allsize += $sizes{$_} if defined $sizes{$_};
  }
  ($bestsize, $best) = choosebest $allsize, join('', @all), $bestsize, $best;
  for my $i ($offset .. $#_) {
    for ($i + 1 .. $#_) {
      $_[$offset] =~ /(\d+)$/;
      next if $_[$_] < $1;
      my @newpartition = @_;
      $newpartition[$i] .= ' ' . splice @newpartition, $_, 1;
      partitionscans $i, @newpartition;
    }
  }
}

my $data = `"$jpegtran" -? 2>&1`;
$tran = $data =~ /outputfile/ ? \&wintran : \&tran;
$data = `"$jpegtran2" -? 2>&1`;
if ($data =~ /outputfile/) {
  $data = `"$jpegtran2" -v @strip -optimize "$fin" "$jtmp" 2>&1`;
  $tran2 = \&wintran;
} else {
  open my $OLDERR, '>&', STDERR;
  open STDERR, '>', $ftmp;
  open my $TRAN, '-|', $jpegtran2, '-v', @strip, '-optimize', $fin;
  $data = <$TRAN>;
  close $TRAN;
  open STDERR, '>&', $OLDERR;
  write_file $jtmp, $data;
  $data = read_file $ftmp;
  $tran2 = \&tran;
}
$data =~ /components=(\d+)/ or die "Couldn't read file";
$planes = $1 - 1;
$width = ($planes + 1) * 2 + 10;
$width2 = length $insize;

for (0 .. $planes) {
  print "Calculating sizes of AC scans, plane $_\n";
  generateacscans $_;
  doscans;
  print "\nSearching for best AC scan, plane $_\n";
  ($acsize[$_], $acscans[$_]) = trysplits $_, 1, 63;
  $acallsize += $acsize[$_];
  $acallscans .= $acscans[$_];
  print "Best AC scan, plane $_:\n$acscans[$_]Size: $acsize[$_]\n\n" if $verbose;
}

print "Calculating sizes and searching for best DC scan\n";
if ($incompat > 0 && $planes > 1) {
  partitionscans 0, 0 .. $planes;
} else {
  partitionscans $planes + 1, 0 .. $planes if $planes > 0 && $incompat >= 0;
  partitionscans 2, join ' ', 0 .. $planes;
}
if ($verbose) {
  my $dcallscans = $best;
  $dcallscans =~ s/$acallscans//;
  print "\nBest DC scan:\n${dcallscans}Size: ", $bestsize - $acallsize, "\n" if $verbose && $dcallscans ne $best;
}

$best = join "\n", sort {
  my $pattern = qr/^([^:]*)(?:: (\d+) (\d+) \d+ (\d+))?/;
  my ($ap, $aa, $ab, $ad) = $a =~ /$pattern/g;
  my ($bp, $ba, $bb, $bd) = $b =~ /$pattern/g;
  !defined $ab || !defined $bb ? $ap cmp $bp : $ab <=> $bb || $ap cmp $bp;
} split "\n", $best;

my @best = map {[/^([^:]*)(?:: (\d+) (\d+) \d+ (\d+))?/g]} split "\n", $best;

my @prev = ();
for my $i (0 .. $#best) {
  for my $j (0 .. $#best) {
    next if $i == $j;
    push @{$prev[$j]}, $i if (
      ($best[$j][0] eq $best[$i][0] && ($best[$j][1] >= $best[$i][1] && $best[$j][1] <= $best[$i][2] || $best[$i][1] >= $best[$j][1] && $best[$i][1] <= $best[$j][2]) && $best[$i][3] > $best[$j][3])
      || ($best[$j][0] gt $best[$i][0] && $best[$j][2] >= $best[$i][2])
    );
  }
}

my $i = -1;
my @bestout = ();
my %bestout = ();
loop: while(@bestout < @best) {
  $i = 0 if ++$i == @best;
  next if defined $bestout{$i};
  for my $j (@{$prev[$i]}) {
    next loop if !defined $bestout{$j};
  }
  push @bestout, $i;
  $bestout{$i} = 1;
  $i = -1;
}
$best = join "\n", (split "\n", $best)[@bestout];

write_file $ftmp, $best;
$data = &$tran('-scans', "\"$ftmp\"", $arith, @strip, "\"$jtmp\"");
$data = app0remove $data if $app0;

if ($jpegtran ne $jpegtran2) {
  $jpegtran = $jpegtran2;
  my $data2 = &$tran2('-scans', "\"$ftmp\"", $arith, @strip, "\"$jtmp\"");
  $data2 = app0remove $data2 if $app0;
  $data = $data2 if length $data2 && length $data2 < length $data;
}

my $size = length $data;
print "\nBest scans:\n$best\nSize: $size bytes\n";
printf "Change: %+.6f%%\n", ($size / $insize - 1) * 100;

if ($size && $size <= $insize) {
  write_file $fout, $data;
} elsif ($fin ne $fout) {
  write_file $fout, read_file $fin;
}
unlink $ftmp, $jtmp;
