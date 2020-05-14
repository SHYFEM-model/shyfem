#!/usr/bin/python

#------------------------------------------------------------------------
#
#    Copyright (C) 1985-2020  Georg Umgiesser
#
#    This file is part of SHYFEM.
#
#------------------------------------------------------------------------

"""Send the contents of a directory as a MIME message.

Usage: dirmail [options] from to [to ...]*

Options:
    -h / --help
        Print this message and exit.

    -d directory
    --directory=directory
        Mail the contents of the specified directory, otherwise use the
        current directory.  Only the regular files in the directory are sent,
        and we don't recurse to subdirectories.

`from' is the email address of the sender of the message.

`to' is the email address of the recipient of the message, and multiple
recipients may be given.

The email is sent by forwarding to your local SMTP server, which then does the
normal delivery process.  Your local machine must be running an SMTP server.
"""


import sys
import os
import getopt
import errno
import mimetypes
import email

from mailbox import UnixMailbox

def print_header(msg,i=0):
    fromh = msg.get("from")
    date = msg.get("date")
    #mydate = '%4d/%2d/%2d' % (date[2],date[1],date[0])
    print i,date,fromh

def scan_file(file):

  fp = open(file, 'r')
  mbox = UnixMailbox(fp)

  #msg = email.message_from_file(fp)
  #print_header(msg)
  #msg = email.message_from_file(fp)
  #print_header(msg)
  #print "ok"
  #return



  i = 0

  for mail in mbox:
    i = i + 1
    fromh = mail.getheader("from")
    date = mail.getdate("date")
    mydate = '%4d/%2d/%2d' % (date[2],date[1],date[0])
    print i,mydate,fromh

def usage(code, msg=''):
    print >> sys.stderr, __doc__
    if msg:
        print >> sys.stderr, msg
    sys.exit(code)


def main():
    try:
        opts, args = getopt.getopt(sys.argv[1:], 'hd:', ['help', 'directory='])
    except getopt.error, msg:
        usage(1, msg)

    dir = os.curdir
    for opt, arg in opts:
        if opt in ('-h', '--help'):
            usage(0)
        elif opt in ('-d', '--directory'):
            dir = arg

    if len(args) < 1:
        usage(1)

    file = args[0]

    scan_file(file)


if __name__ == '__main__':
    main()


