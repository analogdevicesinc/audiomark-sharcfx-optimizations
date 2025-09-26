# http://sourceware.org/gdb/wiki/FAQ: to disable the
# "---Type <return> to continue, or q <return> to quit---"
# in batch mode:
set width 0
set height 0
set verbose off

set logging file output.log
set logging on
show logging


target extended-remote :3333
file ./audiomark-test-abf/audiomark-test-abf.dxe 
load
continue

