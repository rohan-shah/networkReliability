#ifndef EDGE_STATE_HEADER_GUARD
#define EDGE_STATE_HEADER_GUARD
namespace networkReliability
{
	enum EdgeState
	{
		FIXED_OP = 1, FIXED_INOP = 2, UNFIXED_OP = 4, UNFIXED_INOP = 8, NEW_FIXED_INOP = 16, OP_MASK = 5, INOP_MASK = 26, FIXED_MASK = 19, UNFIXED_MASK = 12
	};
}
#endif