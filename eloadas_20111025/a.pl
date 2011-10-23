#!/usr/bin/perl
#----------------------------------------------------------------------------
# "THE BEER-WARE LICENSE" (Revision 42):
# <boroczki.lajos@kfkizrt.hu> wrote this file. As long as you retain this 
# notice you can do whatever you want with this stuff. If we meet some day, 
# and you think this stuff is worth it, you can buy me a beer in return.
# Lajos Boroczki
#----------------------------------------------------------------------------
#
my $pi=3.14167346;

my $facenum=40;
my $levelnum=40;

foreach ((0..$levelnum)){
  my $level=-1+2*$_/$levelnum;
  foreach ((0..$facenum)){
    print sin(2*$_*$pi/$facenum)." ".cos(2*$_*$pi/$facenum)." $level ";
  }
  print "\n";
}
