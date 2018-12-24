#!/bin/bash
#
# creates a release on github given a specific tag
#
# example: git-release.sh VERS_7_5_53 'Release description'
#
# before using this script we need to create a token
#
# https://github.community/t5/How-to-use-Git-and-GitHub/How-to-create-full-release-from-command-line-not-just-a-tag/td-p/6895
#
#---------------------------------------------------------------------

#---------------------------------------------------------------------
# functions
#---------------------------------------------------------------------

Usage()
{
  echo "Usage: git-release.sh [-h|-help] [-options]"
}

FullUsage()
{
  Usage
  echo "Applies release name to exisiting tag"
  echo "Available options:"
  echo "  -h|-help           this help"
  echo "  -write             really write release"
  echo "  -tag tag           use tag to which apply release (default: last)"
  echo "  -title title       use title for the release (default: tag)"
  echo "  -description text  use text as description (default: title)"
  echo "  -draft             release is a draft"
  echo "  -prerelease        release is a prerelease"
  echo "example: git-release.sh VERS_7_5_53 'New Release'" \
		"'This is a new release'"
}

#---------------------------------------------------------------------

CheckVar()
{
  name=$1
  content=$2

  if [ -z "$content" ]; then
    echo "  *** variable $name is empty... error"
    error="YES"
  else
    echo "  $name = $content"
  fi
}

generate_post_data()
{
  cat <<EOF
{
  "tag_name": "$tag",
  "target_commitish": "$branch",
  "name": "$release_name",
  "body": "$release_text",
  "draft": $draft,
  "prerelease": $prerelease
}
EOF
}

#---------------------------------------------------------------------
# parse command line
#---------------------------------------------------------------------

error="NO"
write="NO"
draft="false"
prerelease="false"

while [ -n "$1" ]
do
   case $1 in
        -write)         write="YES";;
        -tag)           tag=$2; shift;;
        -title)         title=$2; shift;;
        -description)   description=$2; shift;;
        -draft)         draft="true";;
        -prerelease)    prerelease="true";;
        -h|-help)       FullUsage; exit 0;;
        -*)             echo "No such option : $1"; exit 1;;
        *)              break;;
   esac
   shift
done
 
if [ $# -ne 0 ]; then
  echo "*** too many arguments on command line: $*"
  Usage
  exit 0
fi

[ -z "$tag" ] && tag=$( make tag )
[ -z "$title" ] && title=$tag
[ -z "$description" ] && description=$title

release_name=$title
release_text=$description

branch=$(git rev-parse --abbrev-ref HEAD)
repo_full_name=$(git config --get remote.origin.url \
		| sed 's/.*:\/\/github.com\///;s/.git$//')
token=$(git config --global github.token)

#---------------------------------------------------------------------
# check available information
#---------------------------------------------------------------------

generate_post_data		# just for check

CheckVar tag $tag
CheckVar release_name "$release_name"
CheckVar release_text "$release_text"
CheckVar branch $branch
CheckVar repo_full_name $repo_full_name
CheckVar token $token
echo

if [ $error = "YES" ]; then
  echo "*** errors found... aborting"
  echo
  Usage
  exit 1
fi

if [ $write = "NO" ]; then
  echo "*** no release written... use -write to really write release"
  echo
  Usage
  exit 0
fi

#---------------------------------------------------------------------
# commit release
#---------------------------------------------------------------------

echo "Create release $tag for repo: $repo_full_name branch: $branch"
curl --data "$(generate_post_data)" \
   "https://api.github.com/repos/$repo_full_name/releases?access_token=$token"

#---------------------------------------------------------------------
# end of routine
#---------------------------------------------------------------------

