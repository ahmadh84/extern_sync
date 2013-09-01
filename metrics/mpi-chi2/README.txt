A short implementation of the chi-squared distance between sets 
of histograms. This is a very useful distance e.g. in the area of 
object recognition.

The code uses compiler intrinsics (on x86 and x86-64) and OpenMP 
parallelism if supported, making it -to my knowledge- the fastest 
implementation that is freely available.

See LICENSE.txt for the licensing terms.

(C) 2007-2008, Christoph Lampert <christoph.lampert@gmail.com>

Changelog:
	2008/07/21: bug fixed in SHUFFLE for double entries.
	Thanks for Sebastian Nowozin for spotting this.

        2009/02/15: bug fixed in symmetric chi2 calculation
                    and python wrapper. Thanks to Nicolas Pinto 
                    for spotting this.
