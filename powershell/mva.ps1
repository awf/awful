
# The classic "mv *.jpg *.jpg.txt" script
# Patterns on the source can be referred to in the destination


# usage function
function usage {
  param([string]$msg)
  write-host @"
usage warning: $msg
usage: mva [-y] [command] src dst

Examples:
   mva *.png *.jpg
   mva *_*.png dir*/file*.png
   mva imconvert *.png *.jpg

Note that without [-y], nothing is moved: 
the script just prints what it *would* do.  
Study that output carefully to make sure it 
does what you expect before re-running with -y
"@
  exit 1
}

if (!$args -or $args.count -eq 0) {
  usage
}

# parse args
$doit = $false;
foreach ($a in $args) {
  if ($a -eq '-y') {
    $doit = $true;
  } else {
    if ($arg3) {
      usage "too many args"
    } elseif ($arg2) {
      $arg3 = $a;
    } elseif ($arg1) {
      $arg2 = $a;
    } else {
      $arg1 = $a
    }
  }
}

if ($arg3) {
  $cmd = $arg1
  $src = $arg2
  $dst = $arg3
} else {
  $cmd = 'mv'
  $src = $arg1
  $dst = $arg2
}

# summarize args
write-debug "Arguments: doit[$doit] cmd[$cmd] src[$src] dst[$dst]"

$srcdir = split-path -parent $src
$srcfile = split-path -leaf $src

if ($srcdir -match '[][*?]') {
  throw "Can't handle globs in the parent directory at the moment"
}

if ($srcdir -eq '') {
  $srcdir = '.';
}

# convert src and dst globs to regexps
$dbp = $DebugPreference
$DebugPreference = 'silentlycontinue'
$srcre = "^" + (awf-glob2regex $srcfile) + "$"
$DebugPreference = $dbp

write-debug "Source glob '$srcfile' becomes regexp '$srcre'"

if ($dst -match '\|') {
  throw "Character '|' not allowed in input"
} 

$dstdir = split-path -parent $dst
$dstfile = split-path -leaf $dst
if ($dstdir -eq '') { $dstdir = '.' }

if ($dstdir -match '[][*?]') {
  throw "Can't handle funny business in the parent directories at the moment"
}

# now replace globs in $dst with backrefs
function replace-backrefs($glob) { 
  $glob = $glob -replace '(\*)|(\?)|({[^}]+})','$|'
  $refindex = 1;
#  write-debug "1: $glob"
  while ($refindex -lt 10) {
    $i = $glob.IndexOf('|');
    if ($i -eq -1) {
      break
    }
    $glob = $glob.SubString(0, $i) + "$refindex" + $glob.SubString($i+1, $glob.length - $i - 1);
    ++$refindex
  }
  $glob
}
# test the above
if (
  ((replace-backrefs 'fr??A{1,2}ed*.*') -ne 'fr$1$2A$3ed$4.$5') -or
  ((replace-backrefs 'pdfs/file***.png') -ne 'pdfs/file$1$2$3.png')
) {
  throw "replace-backrefs failed"
}

$dstpattern = replace-backrefs $dstfile
write-debug "Destination glob '$dstfile' becomes replacement string '$dstpattern'"

# Check for collisions
$d = dir $src
$src_parent_fullname = (get-item $srcdir).PSPath
$d | foreach {
  if ($_.PSParentPath -ne $src_parent_fullname) {
    throw "File [$d] not in the directory [$srcdir]"
  }

  $thissrcfile = $_.name
  if ($thissrcfile -match $srcre) {
    $thisdstfile = [regex]::Replace($thissrcfile, $srcre, $dstpattern);
    $thisdst = $dstdir + '\' + $thisdstfile;
    $thissrc = $srcdir + '\' + $thissrcfile;
    if ($doit) {
      mv $thissrc $thisdst 
    } else {
      write-host "would mv `"${thissrc}`" `"${thisdst}`""
    }
  } else {
    write-host "no match for [$thissrc] to [$srcre]"
  }
}
