#include <boost/date_time/gregorian/gregorian.hpp>

#include <iostream>


//

namespace bg = boost::gregorian;


// ----------------------------------------------------------------------------

//

// ----------------------------------------------------------------------------

int main() {

    // date as simple ISO string

    std::string stringDate("20191110");

    bg::date todayDate(bg::from_undelimited_string(stringDate));


    // 

    bg::date::ymd_type ymdDate = todayDate.year_month_day();

    bg::greg_weekday weekdayDate = todayDate.day_of_week();


    //

    std::cout

        << "Tutorial made with love from Badprog.com :D on "

        << weekdayDate.as_long_string() << " "

        << ymdDate.day << " "

        << ymdDate.month.as_long_string() << " "

        << ymdDate.year

        << std::endl;

}
