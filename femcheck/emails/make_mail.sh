#!/bin/bash
#
#--------------------------------------------------------

mail_file=to_mail.txt

./parse_email.pl addresses.txt > $mail_file

while read -r line; do
    echo "$line"
done < $mail_file 

#--------------------------------------------------------

