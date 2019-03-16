rm *phi1*.root

root -l 'signalRate.cc("zp200phi1")' &>log.zp200phi1 &
root -l 'signalRate.cc("zp300phi1")' &>log.zp300phi1 &
root -l 'signalRate.cc("zp400phi1")' &>log.zp400phi1 &
root -l 'signalRate.cc("zp500phi1")' &>log.zp500phi1 &
root -l 'signalRate.cc("zp750phi1")' &>log.zp750phi1 &
root -l 'signalRate.cc("zp1000phi1")' &>log.zp1000phi1 &
root -l 'signalRate.cc("zp1250phi1")' &>log.zp1250phi1 &
root -l 'signalRate.cc("zp1500phi1")' &>log.zp1500phi1 &
root -l 'signalRate.cc("zp1750phi1")' &>log.zp1750phi1 &
root -l 'signalRate.cc("zp2000phi1")' &>log.zp2000phi1 &

sleep 30;

hadd rate_phi1.root zp*phi1.root

