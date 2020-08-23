function [iborCurr] = getIbor(ibor, iborDates, currDate)

for i = 1:length(iborDates)
   if iborDates(i) == currDate
       ind = i;
       break
   end   
end

iborCurr = ibor(ind);



end