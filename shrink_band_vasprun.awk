#!/usr/bin/awk -f

# Initialise the flags for reading the XML
BEGIN {
	in_dos=0
	in_pdos=0
	in_projected=0
	in_projected_array=0
	in_eigenvalues=0
}
# Read the beginning and end of the XML tags.
{ if ($1 ~ /<projected>/)
	in_projected=1
}
in_projected {
	if ($1 ~ /<\/projected>/)
		in_projected=0
	if ($1 ~ /<eigenvalues>/)
		in_eigenvalues=1
	if (!in_eigenvalues && ($1 ~ /<array>/))
		in_projected_array=1
}
in_eigenvalues {
	if ($1 ~ /<\/eigenvalues>/)
		in_eigenvalues=0
}
in_projected_array {
	if ($1 ~ /<\/array>/)
		in_projected_array=0
}
{ if ($1 ~ /<dos>/)
	in_dos=1
}
in_dos {
	if ($1 ~ /<\/dos>/)
		in_dos=0
	if ($1 ~ /<partial>/)
		in_pdos=1
}
in_pdos {
	if ($1 ~ /<\/partial>/)
		# Unset the flag, and also skip the printing.
		{in_pdos=0; next}
}
# Do special processing of <projected><array>
# Specifically to reduce the precision.
in_projected_array {
	if ($1 ~ /<r>/)
		# Numerical record, do special processing
		printf "        %s %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f %s   \n", $1, $2, $3, $4, $5, $6, $7, $8, $9, $10, $11
	else
		print
}
# Delete the <dos><partial> tag
# Meaning, print when not in it
!(in_pdos) && !(in_projected_array) {print}