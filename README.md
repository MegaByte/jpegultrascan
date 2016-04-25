# jpegultrascan
JPEG recompressor that tries all scan possibilities to minimize size losslessly

## Command line usage
    perl jpegultrascan.pl [switches] inputfile outputfile

## Switches
    -s      Strip all extra markers
    -s some Strip only non-comment markers
    -i      Allow multiple planes per DC scan, which may improve compression
            (incompatible with some software)
    -i -1   Disallow single scan DC components, which may be incompatible with
            some software
    -j      Strip APP0 segment for 18-byte savings (generates non-compliant JPEG
            that may be incompatible with some software)
    -a      Use arithmetic coding (unsupported by most software; jpegtran support
            required)
    -b 0-9  Maximum number of bit splits to test (default: 3)
    -t [N]  Number of simultaneous processes (default: 8 if specified,
            1 otherwise)
    -q      Suppress all output
    -v      Verbose output
    -p path Path to jpegtran (default: jpegtran)
    -o path Path to other jpegtran (default: jpegtran)

## Requirements
- [Perl](https://www.perl.org/)
- [jpegtran](http://www.infai.org/jpeg/) and/or [mozjpeg](https://github.com/mozilla/mozjpeg)

## Notes
jpegultrascan is modeled after [jpegrescan](https://github.com/kud/jpegrescan), but improves compression over jpegrescan by 5-10% by performing an exhaustive search. Resulting optimized images may have many more scans than progressive encoding, though there should be no noticeable difference on modern decoders.

Due to the exhaustive search, recompression takes longer than jpegrescan. Higher values of `-b` results in longer compression times. Few images are likely benefit to from higher `-b` values.

Optimal compression is achieved by specifying the mozjpeg version of jpegtran for the `-o` switch. However, mozjpeg is slower and buggy with certain scan configurations, so IJG jpegtran is recommended for the `-p` switch.

`-i` generates files that are incompatible with some software such as Photoshop and Opera <= 11.61. `-i -1` generates files that are compatible with IPP JPEG.

jpegultrascan is compatible with arbitrary numbers of channels, including CMYK.

## Author
Aaron Kaluszka <<megabyte@kontek.net>>
