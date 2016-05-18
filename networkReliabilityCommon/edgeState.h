#ifndef EDGE_STATE_HEADER_GUARD
#define EDGE_STATE_HEADER_GUARD
namespace networkReliability
{
	enum edgeState
	{
		FIXED_OP = 1, FIXED_INOP = 2, UNFIXED_OP = 4, UNFIXED_INOP = 8, OP_MASK = 5, INOP_MASK = 10, FIXED_MASK = 3, UNFIXED_MASK = 12
	};
}
#endif
