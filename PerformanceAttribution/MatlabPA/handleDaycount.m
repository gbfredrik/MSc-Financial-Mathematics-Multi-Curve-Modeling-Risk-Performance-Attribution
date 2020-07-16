function [dt] = handleDaycount(dcType, daysToCashflow)

if dcType == 'act/360'
    dt = daysToCashflow/360;
end


end