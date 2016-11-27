# jpegdropscan
JPEG lossy recompressor that removes least informative scans to reach a target quality

## Command line usage
    perl jpegdropscan.pl [switches] inputfile outputfile

## Switches
    -d N      Maximum butteraugli difference level (default: 0)
    -r path   Path to reference image to use with butteraugli (default: inputfile)
    -q        Suppress all output
    -v        Verbose output
    -p path   Path to jpegtran (default: jpegtran)
    -c path   Path to butteraugli (default: butteraugli)

## Requirements
- [Perl](https://www.perl.org/)
- jpegtran from [libjpeg](http://www.infai.org/jpeg/), [libjpeg-turbo](https://github.com/libjpeg-turbo/libjpeg-turbo), or [mozjpeg](https://github.com/mozilla/mozjpeg)
- [butteraugli](https://github.com/google/butteraugli)

## Notes
Rather than recoding an image, jpegdropscan iteratively removes scan information starting with least significant bits and moving to least significant coefficients. Butteraugli is used for quality evaluation.

Best results can be achieved by processing an image with [jpegultrascan](jpegultrascan.md) first, running jpegdropscan, then jpegultrascan once more.

## Author
Aaron Kaluszka <<megabyte@kontek.net>>
