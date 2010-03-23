#!/bin/sh

# gets headers of fortran files

echo "<html><head><title>Fortran Headers</title></head><body>"

echo ""
echo "<h1>"
echo "Fortran Headers"
echo "</h1>"
echo ""
echo ""

files=$*

echo "<h2><a name=\"index\">Index</a></h2>"

for file in $files
do
  echo "<a href=\"#$file\">$file</a>"
done

echo "<h2>"
echo "Headers"
echo "</h2>"
echo ""
echo ""

for file in $files
do
  echo "<h3><a name=\"$file\">"
  echo "$file"
  echo "</a></h3>"
  echo ""
  echo "<pre>"
  getheader.pl $file
  echo "</pre>"
  echo "<small>Back to <a href=\"#index\">index</a></small>"
  echo ""
done

echo "</body></html>"
