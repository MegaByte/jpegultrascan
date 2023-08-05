#!/usr/bin/perl
=head1 NAME
jpegdropscan
=head1 DESCRIPTION
JPEG lossy recompressor that removes least informative scans to reach a target quality
=head1 VERSION
1.0.2 2023-06-10
=head1 LICENSE
Copyright 2015 - 2023 Aaron Kaluszka

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

my $jpegtran = 'jpegtran';
my $butteraugli = 'butteraugli';
my $home = glob '~';
my $quality = 0;
my $maxbits = 3;
GetOptions 'q' => \my $quiet, 'v' => \my $verbose, 'd=f' => \$quality,
           'p=s' => \$jpegtran, 'c=s' => \$butteraugli, 'r=s' => \my $reference, 'h|?' => \my $help;
$verbose = defined $verbose;
$jpegtran =~s /^~/$home/;
$butteraugli =~s /^~/$home/;
$quality = 0 if $quality < 0;

open STDOUT, '>', File::Spec->devnull if defined $quiet;

@ARGV == 2 && !defined $help or die 'usage: jpegdropscan.pl [switches] inputfile outputfile
JPEG lossy recompressor that removes least informative scans to reach a
target quality
Switches:
  -q        Suppress all output
  -v        Verbose output
  -p path   Path to jpegtran (default: jpegtran)
  -c path   Path to butteraugli (default: butteraugli)
  -r path   Path to reference image to use with butteraugli (default: inputfile)
  -d N      Maximum butteraugli difference level (default: 0)
';

my ($fin, $fout) = @ARGV;
$reference = $fin unless defined $reference;
$fin =~s /^~/$home/;
$fout =~s /^~/$home/;
$reference =~s /^~/$home/;
my ($tran, @best, $bestsize, $bestscore, $app0);
my %colormap = (
  1 => 0, 2 => 1, 3 => 2, 4 => 3, # YCC(K)
  82 => 0, 71 => 1, 66 => 2, # RGB
  67 => 0, 77 => 1, 89 => 2, 75 => 3, # CMYK
  114 => 0, 103 => 1, 98 => 2, # BG RGB
  34 => 1, 35 => 2 # BG YCC
);
my $mapcolor = \&mapcolor;
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
  read_file $tmp;
}

sub butteraugli($) {
  chop(my $score = `"$butteraugli" "$reference" "$_[0]"`);
  $score;
}

sub mapcolor2($) {
  $_[0];
}

sub mapcolor($) {
  if ($_[0] == 0) {
    $mapcolor = \&mapcolor2;
    return 0;
  }
  $colormap{$_[0]};
}

sub scaninfo($) {
  my $file = substr $_[0], 2;
  $file = substr $file, $+[0] while $file =~ /\xFF\xD8.+?\xFF\xD9/s; # skip thumbnail(s)
  my @scans;
  while ($file =~ /\xFF\xDA/gs) {
    my $ns = ord substr $file, $+[0] + 2, 1;
    my $cs = join ' ', map &$mapcolor($_), unpack "(Cx)$ns", substr $file, $+[0] + 3, $ns * 2;
    my ($ss, $se, $ahal) = unpack 'C3', substr $file, $+[0] + 3 + $ns * 2, 3;
    push @scans, [$cs, $ss, $se, $ahal >> 4, $ahal & 3];
  }
  @scans;
}

sub app0remove($) {
  my ($file) = @_;
  while ("\xFF\xE0" eq substr $file, 2, 2) {
    my $size = 2 + unpack 'n', substr $file, 4, 2;
    substr $file, 2, $size, '';
  }
  $file;
}

sub printscans(@) {
  join '', map sprintf("%s: %d %d %d %d;\n", @$_), @_;
}

sub tryscans($$) {
  my ($file, $scans) = @_;
  my (undef, $stmp) = tempfile SUFFIX => $$;
  write_file $stmp, $scans;
  my $data = &$tran('-scans', "\"$stmp\"", "\"$file\"");
  $data = app0remove $data if $app0;
  $data;
}

sub clone(@) {
  map [@$_], @_;
}

sub checkbest(@) {
  my $scans = printscans @_;
  print "\n$scans" if $verbose;
  my (undef, $tmp) = tempfile SUFFIX => $$;
  my $result = tryscans $fin, $scans;
  write_file $tmp, $result;
  my $score = butteraugli $tmp;
  print $verbose ? "$score: " . length($result) . "\n" : '.';
  if ($score > $quality) {
    return -1;
  } elsif (length $result < $bestsize || length $result == $bestsize && $score >= $bestscore && @_ <= @best) {
    @best = clone @_;
    $bestsize = length $result;
    $bestscore = $score;
    return 1;
  }
  0;
}

sub inside($$$$) {
  $_[0] >= $_[2] && $_[0] <= $_[3] || $_[2] >= $_[0] && $_[2] <= $_[1];
}

sub shavebits() {
  my @s = clone(@best);
  for (my $b = 1; $b <= $maxbits; ++$b) {
    attempt: for (my $s = $#s; $s > 0; --$s) {
      my @lastscan = clone @s;
      next if $s[$s][4] != $b - 1;
      if ($s[$s][3] != 0) {
        for (my $t = $s + 1; $t < @s; ++$t) {
          next attempt if $s[$s][0] eq $s[$t][0] && inside($s[$s][1], $s[$s][2], $s[$t][1], $s[$t][2]) && $s[$t][3] < $b;
        }
        splice @s, $s, 1;
        @s = @lastscan if checkbest(@s) == -1;
      } else {
        $s[$s][4] = $b;
        my $skip = 0;
        for (my $t = $s + 1; $t < @s; ++$t) {
          if ($s[$s][0] eq $s[$t][0] && inside($s[$s][1], $s[$s][2], $s[$t][1], $s[$t][2]) && $s[$t][3] < $b) {
            $skip = 1;
            last;
          }
        }
        if ($skip || checkbest(@s) == -1) {
          my ($p, $q) = @{$s[$s]}[1, 2];
          if ($p < $q) {
            $s[$s][4] = $b - 1;
            splice @s, $s + 1, 0, [$s[$s][0], $p, $q, $s[$s][3], $b];
            subattempt: for (; $p < $q; ++$p) {
              $s[$s][2] = $p;
              $s[$s + 1][1] = $p + 1;
              for (my $t = $s + 2; $t < @s; ++$t) {
                next subattempt if $s[$s][0] eq $s[$t][0] && $s[$t][3] < $b && (inside($s[$s + 1][1], $s[$s + 1][2], $s[$t][1], $s[$t][2]) || inside($s[$s][1], $s[$s][2], $s[$t][1], $s[$t][2]));
              }
              next attempt if checkbest(@s) == 1;
            }
          }
          @s = @lastscan;
        }
      }
    }
  }
}

sub shavecoefs() {
  my @s = clone @best;
  for (my $c = 62; $c >= 0; --$c) {
    attempt: for (my $s = $#s; $s > 0; --$s) {
      next if $s[$s][2] != $c + 1;
      for (my $t = $s + 1; $t < @s; ++$t) {
        next attempt if $s[$s][0] eq $s[$t][0] && inside($s[$s][1], $s[$s][2], $s[$t][1], $s[$t][2]) && $s[$t][2] > $c;
      }
      my @lastscan = clone @s;
      if ($s[$s][1] > $c) {
        splice @s, $s, 1;
        @s = @lastscan if checkbest(@s) == -1;
      } elsif ($s[$s][2] > $c) {
        $s[$s][2] = $c;
        @s = @lastscan if checkbest(@s) == -1;
      }
    }
  }
}

my $origdata = read_file $fin;
die "$fin is not a valid JPEG file" if length($origdata) < 4 || "\xFF\xD8" ne substr $origdata, 0, 2;
$app0 = "\xFF\xE0" eq substr $origdata, 2, 2;
$tran = `"$jpegtran" -? 2>&1` =~ /outputfile/ ? \&wintran : \&tran;
my @scans = sort {
  my ($ap, $aa, $ab) = @$a;
  my ($bp, $ba, $bb) = @$b;
  !defined $ab || !defined $bb ? $ap cmp $bp : $ab <=> $bb || $ap cmp $bp;
} scaninfo $origdata;

my @prev = ();
for my $i (0 .. $#scans) {
  for my $j (0 .. $#scans) {
    next if $i == $j;
    push @{$prev[$j]}, $i if ($scans[$j][0] eq $scans[$i][0] && inside($scans[$i][1], $scans[$i][2], $scans[$j][1], $scans[$j][2]) && $scans[$i][4] > $scans[$j][4]) || ($scans[$j][0] gt $scans[$i][0] && $scans[$j][2] >= $scans[$i][2]);
  }
}

my $i = -1;
my @order = ();
my %order = ();
loop: while(@order < @scans) {
  $i = 0 if ++$i == @scans;
  next if defined $order{$i};
  for my $j (@{$prev[$i]}) {
    next loop unless defined $order{$j};
  }
  push @order, $i;
  $order{$i} = 1;
  $i = -1;
}
@scans = @scans[@order];

for (my $i = 0; $i < @scans; ++$i) {
  my @s = @{$scans[$i]};
  if ($s[1] == 0 && $s[2] > 0) {
    splice @scans, $i++, 1, [$s[0], 0, 0, $s[3], $s[4]];
    splice @scans, $i++, 0, [$_, 1, $s[2], $s[3], $s[4]] for split ' ', $s[0];
  }
}
my $origscore = butteraugli $fin;
my $insize = -s $fin;
print printscans(@scans), "$origscore: $insize\n" if $verbose;
die "Initial quality $origscore is less than desired quality." if $origscore > $quality;

$bestsize = length $origdata;
$bestscore = $origscore;
@best = @scans;
shavebits;
shavecoefs;

my $best = printscans @best;
my $data = tryscans $fin, $best;
print "\nBest scans:\n$best\nSize: $bestsize bytes\n";
printf "Change: %+.6f%%\n", ($bestsize / $insize - 1) * 100;
print "Quality: $bestscore\n";
printf "Change: %+.6f\n", ($bestscore - $origscore);
printf "Total Change: %+.6f%%\n", ($bestsize / (-s $reference) - 1) * 100 if $reference ne $fin;

if ($bestsize && $bestsize <= $insize) {
  write_file $fout, $data;
} elsif ($fin ne $fout) {
  write_file $fout, $origdata;
}
