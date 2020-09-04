/**
 * @brief Accumulator structure for the Hough transform
 * @file Accumulator.h
 *
 * This file is created at Almende B.V. It is open-source software and part of the Common
 * Hybrid Agent Platform (CHAP). A toolbox with a lot of open-source tools, ranging from
 * thread pools and TCP/IP components to control architectures and learning algorithms.
 * This software is published under the GNU Lesser General Public license (LGPL).
 *
 * It is not possible to add usage restrictions to an open-source license. Nevertheless,
 * we personally strongly object to this software being used by the military, in factory
 * farming, for animal experimentation, or anything that violates the Universal
 * Declaration of Human Rights.
 *
 * Copyright © 2013 Anne van Rossum <anne@almende.com>
 *
 * @author  Anne C. van Rossum
 * @date    May 7, 2013
 * @project Replicator FP7
 * @company Almende B.V.
 * @case    Artificial Intelligence Framework
 */

#ifndef ACCUMULATOR_H_
#define ACCUMULATOR_H_

#include <HoughDefs.h>

#include <nd-array.hpp>

namespace dobots {

/**
 * The Hough transform uses a specific type of data structure that is called the accumulator in this context. It is
 * basically just a discretization of the Hough space. For lines, which have only two parameters, you will need a 2D
 * Hough space.
 */
template <typename P>
class Accumulator : public nd_array<Cell<P>,ACCUMULATOR_DATA_TYPE> {
public:
	//! Default accumulator constructor
	Accumulator(ASize size) {
		this->asize.x = size.x;
		this->asize.y = size.y;
		std::vector<ACCUMULATOR_DATA_TYPE> dimensions; dimensions.clear();
		dimensions.push_back(size.x);
		dimensions.push_back(size.y);
		this->init(dimensions);
	}

	//! Default accumulator destructor
	~Accumulator() {}

	/**
	 * Increment cell in the accumulator. This assumes that the calling function knows how to map from the point cloud
	 * to the coordinates of a given cell. The original points are also handed over to keep track of segments.
	 */
	void Increment(ACoordinates c, Segment2D<P> seg) {
		assert (c.x >= 0 && c.x < asize.x);
		assert (c.y >= 0 && c.y < asize.y);

		Cell<P> & cell = this->get(c.x, c.y);
		cell.hits++;
		cell.segments.push_back(seg);
		cell.points.push_back(seg.src);
		cell.points.push_back(seg.dest);
	}

	/**
	 * Reset the accumulator, for example when switching from method, but also in some methods you are required to do
	 * so.
	 */
	void Reset() {}

	//! Size of the accumulator
	inline ASize getSize() { return asize; }
private:

	ASize asize;
};

}

#endif /* ACCUMULATOR_H_ */
