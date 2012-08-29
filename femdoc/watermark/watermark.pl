#!/usr/bin/perl

$flag = 0;
while (<>) {
    if (/^%%Page:/) {
        if ($flag) {
            print "grestore\n";
        }
        $flag = 1;
        print $_;
        print "gsave\n";
        print ".95 setgray\n";
        print "/Helvetica-Bold findfont 72 scalefont setfont\n";
        print "80 80 800 { 306 exch moveto\n";
        print "(Draft) dup\n";
        print "stringwidth pop 2 div neg 0 rmoveto show } for\n";
        print "grestore\n";
    } else {
        print;
    }
}
