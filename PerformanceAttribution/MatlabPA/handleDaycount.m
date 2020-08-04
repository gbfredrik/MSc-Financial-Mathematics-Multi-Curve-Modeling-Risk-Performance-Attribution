function [dt] = handleDaycount(dcType, daysToCashflow)

if dcType == 'act/360'
    dt = daysToCashflow/365;
end


end