#ifndef COLORIZE_HPP
#define COLORIZE_HPP

#include<termcolor/termcolor.hpp>

template<typename realType> void colorizedInput(realType& param, std::string msg)
{
	std::cout<<termcolor::bright_yellow<<termcolor::on_grey;
	std::cout<<"Enter a value for "<<msg<<termcolor::reset;
	
	std::cout<<termcolor::bright_cyan<<termcolor::on_grey;
	std::cin>>param;
	std::cout<<termcolor::reset;

	std::cout<<termcolor::bright_cyan<<termcolor::on_grey;
	std::cout<<"[RECEIVED] "<<msg<<param<<termcolor::reset<<std::endl;
	std::cout<<std::endl;
}

void colorizedHeading(std::string msg)
{
	std::cout<<termcolor::bright_yellow<<termcolor::on_grey;
	std::cout<<msg<<std::endl;
	std::cout<<termcolor::bright_magenta<<termcolor::on_grey;
	for(int i=0; i<msg.size(); i++)
	{
		std::cout<<"-";
	}
	std::cout<<termcolor::reset<<std::endl;
}

template<typename realType> void colorizedNumericOutput(realType value)
{
	std::cout<<termcolor::bright_cyan<<termcolor::on_grey;
	std::cout<<value<<termcolor::reset<<std::endl;
}

void colorizedStringOutput(std::string msg)
{
	std::cout<<termcolor::bright_yellow<<termcolor::on_grey;
	std::cout<<msg<<termcolor::reset;
}
#endif
