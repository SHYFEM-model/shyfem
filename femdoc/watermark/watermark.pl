#!/usr/bin/perl

$flag = 0;
$text = "Draft";
$text = "Draft ISMAR-CNR";
$text = "ISMAR-CNR";
$gray = 0.95;
$gray = 0.90;
$ymax = 800;
$ymax = 750;
$dy = 80;
$dy = 160;
$y0 = 80;

while (<>) {
    if (/^%%Page:/) {
        if ($flag) {
            print "grestore\n";
        }
        $flag = 1;
        print $_;
        print "gsave\n";
        print "$gray setgray\n";
        print "/Helvetica-Bold findfont 72 scalefont setfont\n";
        print "$y0 $dy $ymax { 306 exch moveto\n";
        print "($text) dup\n";
        print "stringwidth pop 2 div neg 0 rmoveto show } for\n";
        print "grestore\n";
    } else {
        print;
    }
}
