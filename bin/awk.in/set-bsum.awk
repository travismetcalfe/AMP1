#  awk script to extract salient information from csum
BEGIN{n=0}
/^#/ {next}
n == 0 {print "# Input file: ", FILENAME; print "#";
	print "# n, M/Msun, age (Gyr), R, Teff, L/Lsun, Xc, qc"; print "#";
	format="%5d %7.3f %8.5f %13.5e %7.1f %8.3f %12.5e %11.4e\n";
	lsun=3.846e33}
{n += 1; printf format,n, $1, $2/1.e9, $3, $4, $5/lsun, $9, $16}
