#ifndef LOGGER_H
#define LOGGER_H

#include <iostream>
/**
* Set this to <code>true</code> if debug information should be printed during sampling\n
* Logging is by default directed to the standard output but its quite easy to modify this code and
* redirect the output.
*/
constexpr bool DEBUG = false;
/**
* Set this to <code>true</code> if integrity of data structures should be tested during sampling
*/
constexpr bool TEST = false;

/**
* Set this to <code>true</code> if debug information should be printed during sampling\n
* Logging is by default directed to the standard output but its quite easy to modify this code and
* redirect the output.
*/
constexpr bool DEBUG_PARAMETERS = false;

void print();

void printErr();

template <class A0, class ...Args> void print(A0 a0, Args ...args)
{
	std::cout << a0;
	print(args...);
}

template <class A0, class ...Args> void printErr(A0 a0, Args ...args)
{
	std::cerr << a0;
	printErr(args...);
}

template <class ...Args> void log(Args ...args)
{
	return print(args...);
}

template <class TestClass, class ...Args> void logTest(Args ...args)
{
	return print("Test: ", typeid(TestClass).name(), "\n    Message: ", args...);
}

template <class ...Args> void logErr(Args ...args)
{
	return printErr(args...);
}

inline void print()
{
	std::cout << "\n";
}

inline void printErr()
{
	std::cerr << "\n";
}

#endif // !LOGGER_H
