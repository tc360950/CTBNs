#ifndef STATISTICS_H
#define STATISTICS_H
template <class Real_t> struct Statistics {
	size_t d;
	Real_t t_max;
	Real_t power;
	Real_t FDR;
	Real_t MD;
};
#endif // !STATISTICS_H
