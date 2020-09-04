/**
 * @file DetectLineModuleExt.cpp
 * @brief DetectLineModule extension
 *
 * This file is created at Almende B.V. It is open-source software and part of the Common 
 * Hybrid Agent Platform (CHAP). A toolbox with a lot of open-source tools, ranging from 
 * thread pools and TCP/IP components to control architectures and learning algorithms. 
 * This software is published under the GNU Lesser General Public license (LGPL).
 *
 * It is not possible to add usage restrictions to an open-source license. Nevertheless,
 * we personally strongly object against this software being used by the military, in the
 * bio-industry, for animal experimentation, or anything that violates the Universal
 * Declaration of Human Rights.
 *
 * Copyright © 2013 Your Name <your@email>
 *
 * @author	Your Name
 * @date	Current date
 * @company	Your Company
 * @project	Specific Software Project
 */

#include <algorithm>
#include <functional>

#include <DetectLineModuleExt.h>
#include <Random.h>

#include <HoughDefs.h>
#include <cmath>
#include <fstream>

using namespace rur;
using namespace dobots;

//! Adjust these figures, make sure patch width times number of patches unequivocally matches the entire width
const int NR_PATCHES_WIDTH = 16;
const int NR_PATCHES_HEIGHT = 16;

const int TOTAL_TICK_COUNT = 10000;

// threshold is now at at least 1 percent of the total number of points
const int HIT_THRESHOLD = TOTAL_TICK_COUNT / 100;

const float SKEW_THRESHOLD = 0.9; // requires at least 10 points
const int MATCH_NEIGHBOURHOOD_DEVIATION = 15;

//static int DEBUG_DECPOINT_INDEX = 0;

DecPoint::DecPoint(): Point2D(), skew(0), match(NULL) {
}

DecPoint::DecPoint(int x, int y): Point2D(x,y), skew(0), match(NULL) {
}

DetectLineModuleExt::DetectLineModuleExt(): tick(0), stop(false), image(NULL),
		segmentation(), data_item(), hough(NULL) {

//	data_item = DI_PENTAGON;
	data_item = DI_SQUARE;
//	data_item = DI_BOOKS;
	data_item = DI_POINT_LINES;

	segmentation = GLUE_POINTS;
	segmentation = ALL_POINTS;
//	segmentation = LONGEST_LINE;
//
	srand(time(NULL));
}

DetectLineModuleExt::~DetectLineModuleExt() {

}

/**
 * Helper function that considers all points above zero in an image to be points detected previously by a corner or
 * feature detector. If all points are thrown in a vector, it is hard to quickly get out points that are in each
 * other's proximity.
 */
void DetectLineModuleExt::loadImage(std::string file, std::vector<Point2D> & points) {
	image = new CRawImage(640,480,1);
	if (!image->loadBmp(file.c_str())) {
		std::cerr << "Failed to open file " << file << std::endl;
		return;
	}
	if (!image->isMonochrome()) {
		std::cout << "Make monochrome" << std::endl;
		image->makeMonochrome();
	}
	image->saveBmp("test.bmp");

	// add all points to linear array of points, skipping the margins which are most often just artifacts
	int margin = 2;
	for (int i = margin; i < (image->getwidth() - margin); ++i) {
		for (int j = margin; j < (image->getheight() - margin); ++j) {
			if (image->getValue(i,j) > 0) {
				points.push_back(Point2D(i,j));
			}
		}
	}
	std::cout << "Got " << points.size() << " points out of image " << file << std::endl;
}

/**
 * Assumes that the dimensions of the nd-array are previously set. For now we only expect it to be 2D similar to the
 * dimension of the physical pixel array. It becomes "circular" to define proximity in a way that already "knows" their
 * physical 3D location.
 */
void DetectLineModuleExt::loadImage(std::string file, pointcloud & spatial_points) {
	image = new CRawImage(640,480,1);
	if (!image->loadBmp(file.c_str())) {
		std::cerr << "Failed to open file " << file << std::endl;
		return;
	}
	if (!image->isMonochrome()) {
		std::cout << "Make monochrome" << std::endl;
		image->makeMonochrome();
	}
	image->saveBmp("copyofimage.bmp");

	assert (spatial_points.get_dimensions() == 2);

	default_patch.width = image->getwidth() / spatial_points.get_dimension(0);
	default_patch.height = image->getheight() / spatial_points.get_dimension(1);
	std::cout << "Patch size is " << default_patch.width << "x" << default_patch.height << std::endl;

	std::cout << "Make sure the image can be neatly divided into patches, or adjust NR_PATCHES_*" << std::endl;
	assert(default_patch.width * spatial_points.get_dimension(0) == image->getwidth());
	assert(default_patch.height * spatial_points.get_dimension(1) == image->getheight());

	// Loop over all patches and add points to the patch where they belong (this to accelerate search later on)
	for (int x = 0; x != NR_PATCHES_WIDTH; ++x) {
		for (int y = 0; y != NR_PATCHES_HEIGHT; ++y) {
			// Get patch at (x,y)
			std::vector<DecPoint*> & patch = spatial_points.get(x,y);
			patch.clear();
			for (int i = 0; i < default_patch.width; ++i) {
				int xi = x * default_patch.width + i;
				for (int j = 0; j < default_patch.height; ++j) {
					int yj = y * default_patch.height + j;
					// Get value at (xi,yj). If the value is non-zero add it to this patch.
					if (image->getValue(xi,yj) > 0) {
						patch.push_back(new DecPoint(xi,yj));
					}
				}
			}
		}
	}

	// For debugging purposes, count the total number of points.
	long int pnt_count = 0;
	for (int x = 0; x != NR_PATCHES_WIDTH; ++x) {
		for (int y = 0; y != NR_PATCHES_HEIGHT; ++y) {
			std::vector<DecPoint*> & patch = spatial_points.get(x,y);
			pnt_count += patch.size();
		}
	}
	std::cout << "Got " << pnt_count << " points out of image " << file << std::endl;
}

/**
 *
 */
void DetectLineModuleExt::Init(std::string & name) {
#ifdef LINEAR
	std::vector<DecPoint> points;
	loadImage("../data/dots.bmp", points);
	hough.addPoints(points);
#endif
	// create an n-dim array (currently used 2D)
	std::vector<short> dimensions;
	dimensions.push_back(NR_PATCHES_WIDTH);
	dimensions.push_back(NR_PATCHES_HEIGHT);
	pointcloud cloud(dimensions);
	switch (data_item) {
	case DI_SQUARE: default:
		loadImage("../data/square.bmp", cloud);
		break;
	case DI_PENTAGON:
//		loadImage("../data/pentagon.bmp", cloud);
		loadImage("../data/lena.bmp", cloud);
		break;
	case DI_BOOKS:
		loadImage("../data/books.bmp", cloud);
		break;
	case DI_POINT_LINES:
		loadPoints(name, cloud);
		//loadPoints("../data/dpmline4.pnts.Sat_Aug_29_19:28:24_2015.811392.data.txt", cloud);
		//writeImage("../data/dpmline4.pnts.Sat_Aug_29_19:28:24_2015.811392.data.txt", "../data/output.bmp");
		break;
	}
	ISize input_size (image->getwidth(), image->getheight());
	hough = new Hough<DecPoint>(input_size);

	hough->addPoints(cloud);

	tick = 0;
	std::cout << "Loaded all points (normally from image processing module)" << std::endl;
}

void DetectLineModuleExt::writeImage(std::string ifilename, std::string ofile) {

}

void DetectLineModuleExt::loadPoints(std::string ifilename, pointcloud & cloud) {
	struct IndexedPoint {
		double x;
		double y;
		int index;
		int cluster_index;
	};
	std::vector<IndexedPoint> points;

	std::ifstream ifs(ifilename);
	int i = 0;
	if (ifs.is_open()) {
		while(true) {
			double x, y, cluster_index;
			ifs >> x >> y >> cluster_index;
			if (!ifs.good()) break;
			IndexedPoint pnt;
			pnt.x = x;
			pnt.y = y;
			pnt.index = i++;
			pnt.cluster_index = cluster_index;
			points.push_back(pnt);
		}
	} else {
		std::cerr << "File cannot be opened: " << ifilename << std::endl;
	}
	ifs.close();

	std::cout << "Found " << points.size() << " points" << std::endl;

	double min_x = points[0].x, min_y = points[0].y;
	double max_x = min_x, max_y = min_y;
	for (int i = 0; i < points.size(); ++i) {
		//std::cout << points[i].x << '\t' << points[i].y << std::endl;
		if (points[i].x < min_x) min_x = points[i].x;
		if (points[i].x > max_x) max_x = points[i].x;
		if (points[i].y < min_y) min_y = points[i].y;
		if (points[i].y > max_y) max_y = points[i].y;
	}

	// scale x coordinate to 640 and y coordinate to 480
	double dx = max_x - min_x;
	double dy = max_y - min_y;
	for (int i = 0; i < points.size(); ++i) {
		points[i].x = (int)((points[i].x - min_x)*639/dx);
		points[i].y = (int)((points[i].y - min_y)*479/dy);
		/*
		int x = points[i].x;
		int y = points[i].y;
		if (x < 0) exit(1);
		if (y < 0) exit(1);
		if (x > 640) exit(1);
		if (y >= 480) exit(1);
		*/
		//std::cout << points[i].x << '\t' << points[i].y << std::endl;
	}

	int w = 640 / NR_PATCHES_WIDTH;
	int h = 480 / NR_PATCHES_HEIGHT;
	std::cout << "Width: " << w << std::endl;
	std::cout << "Height: " << h << std::endl;
	
	default_patch.width = w;
	default_patch.height = h;
	
	// store in pointcloud
	for (int i = 0; i < points.size(); ++i) {
		//int x = points[i].x / w;
		//int y = points[i].y / h;
		int p_x = points[i].x / w;
		int p_y = points[i].y / h;
		int x = points[i].x;
		int y = points[i].y;
		//std::cout << points[i].x << " -> " << x << '\t' << points[i].y << " -> " << y << std::endl;
		std::vector<DecPoint*> & patch = cloud.get(p_x,p_y);
		DecPoint * pnt = new DecPoint(x,y);
		pnt->index = points[i].index;
		pnt->cluster_index = points[i].cluster_index;
		patch.push_back(pnt);
	}
	std::cout << "Filled point cloud" << std::endl;
	
	// make fake image
	image = new CRawImage(640,480,1);
	if (!image->isMonochrome()) {
		std::cout << "Make monochrome" << std::endl;
		image->makeMonochrome();
	}
	for (int i = 0; i < points.size(); ++i) {
		int x = points[i].x;
		int y = points[i].y;
		image->setValue(x,y,255);
	}
	image->saveBmp("test2.bmp");
}

void DetectLineModuleExt::Tick() {
	//std::cout << "Hough transform tick" << std::endl;
	hough->doTransform();

	if (++tick == TOTAL_TICK_COUNT) {

		std::cout << "Plot the accumulator data structure" << std::endl;
		plotAccumulator();
		
		std::cout << "Write assignments" << std::endl;
		writeAssignments("output.txt");

		std::cout << "Use segmentation type: " << SegmentationDescription[segmentation] << std::endl;

		std::cout << "Prepare the information to be able to do the segmentation" << std::endl;
		prepareSegments();

		std::cout << "Get the actual segments" << std::endl;
		getSegments();

		std::cout << "Plot the segments which should more or less reconstruct the original image" << std::endl;
		plotSegments();

		stop = true;
	}
}

/**
 * Reasoning on the points that contributed to a cell in the accumulator. The segments should be added to the vector
 * given as argument.
 *
 * So here every segment, that is every pair of points, that increments the accumulator at that cell, will be added 
 * to this cell. A segment in e.g. ALL_POINTS is nothing more than each pair of points! If you draw the result, it 
 * can still look good, because all those segments are overlapping on a line.
 */
void DetectLineModuleExt::addSegments(Cell<DecPoint> & c, std::vector<Segment2D<DecPoint> > & sgmnts) {
	if (c.segments.empty()) return;
	switch(segmentation) {
	case ALL_POINTS: default: {
		// just consider every segment how it happens to enter as the final segment, no calculations here, neither
		// interpolations
		sgmnts.insert(sgmnts.end(), c.segments.begin(), c.segments.end());
		break;
	}
	case LONGEST_LINE: {
		//sort(c.points.begin(), c.points.end(), less_than_x());
		// returns only one line segment, the one that is largest, support for this line by the points is not respected
		DecPoint *min_x = c.points[0], *min_y = c.points[0], *max_x = c.points[0], *max_y = c.points[0];
		for (int i = 0; i < c.points.size(); ++i) {
			if (c.points[i]->x < min_x->x) min_x = c.points[i];
			if (c.points[i]->x > max_x->x) max_x = c.points[i];
			if (c.points[i]->y < min_y->y) min_y = c.points[i];
			if (c.points[i]->y > max_y->y) max_y = c.points[i];
		}
		int dx = max_x->x - min_x->x;
		int dy = max_y->y - min_y->y;
		Segment2D<DecPoint> segment;
		if (dx > dy) {
			segment.src = min_x;
			segment.dest = max_x;
		} else {
			segment.src = min_y;
			segment.dest = max_y;
		}
		sgmnts.push_back(segment);
		break;
	}
	case DISTRIBUTION_SENSITIVE: {
		// suppose we have points in the following pattern:
		//   .... .... ... . . .                 ... .... . . .... . . . . . .. .....
		// then we likely want to get two segments out of here (there is no ground truth here)
		// this corresponds to k-clustering or even nonparametric (no k parameter) clustering

		break;
	}
	case GLUE_POINTS: {
		// suppose we have points in the following pattern:
		//   .... .... ... . . .                 ... .... . . .... . . . . . .. .....
		// we would like to obtain the extremities, points that are not "surrounded" by other points

		Segment2D<DecPoint> segment;
		for (int i = 0; i < c.points.size(); ++i) {
			// extreme form, only consider points that have a matching "partner", we do not care about this match,
			// but use it only to indicate that this point is probably the end of a segment...
			if (c.points[i]->match != NULL) {
				if (segment.src == NULL)
					segment.src = c.points[i];
				else {
					if (segment.src == c.points[i]) continue;
					if ((std::abs(segment.src->x - c.points[i]->x) < MATCH_NEIGHBOURHOOD_DEVIATION) &&
							(std::abs(segment.src->y - c.points[i]->y) < MATCH_NEIGHBOURHOOD_DEVIATION))
						continue;
					segment.dest = c.points[i];
					std::cout << "Found segment " << *(segment.src) << " - " << *(segment.dest) << std::endl;
					sgmnts.push_back(segment);
					break;
				}
			}
		}

		// return those points that have a high value for EOS
		break;
	}
	} // end of switch statement
}

/**
 * Get segments by looping over each cell and calling getSegments for each of them if they are beyond a certain
 * threshold w.r.t. number of hits.
 */
void DetectLineModuleExt::getSegments() {
	Accumulator<DecPoint> & acc = *hough->getAccumulator();
	ASize asize = acc.getSize();

	// loop over the accumulator
	for (int i = 0; i < asize.x; ++i) {
		for (int j = 0; j < asize.y; ++j) {
			Cell<DecPoint> & c = acc.get(i,j);
			if (c.hits > HIT_THRESHOLD) {
				addSegments(c, segments);
			}
		}
	}
}

void DetectLineModuleExt::calculateSkew(Cell<DecPoint> * c) {
	if (c->hits <= HIT_THRESHOLD) return;

	sort(c->points.begin(), c->points.end(), xref_increasing());

	for (int i = 0; i < c->points.size(); ++i) {
		float midway = (float) c->points.size() / 2;
		c->points[i]->skew = std::fabs(i - midway) / midway;
		//		std::cout << "Calculated skew for point " << i << " is " << c->points[i]->skew << std::endl;
	}
	sort(c->points.begin(), c->points.end(), skewref_decreasing());
}

void DetectLineModuleExt::findMatches(Cell<DecPoint> * self, Cell<DecPoint> * other) {
	if (self->hits <= HIT_THRESHOLD) return;

	for (int i = 0; i < self->points.size(); ++i) {
		if (self->points[i]->skew < SKEW_THRESHOLD) break; // break, because ordered on skew
		//		std::cout << "Skew first point is above threshold: " << self->points[i]->skew << std::endl;
		for (int j = 0; j < other->points.size(); ++j) {
			if (other->hits <= HIT_THRESHOLD) continue;
			if (other->points[j]->skew < SKEW_THRESHOLD) break; // break, because ordered on skew
			//			std::cout << "Found potential match!" << std::endl;
			if ((std::abs(self->points[i]->x - other->points[j]->x) < MATCH_NEIGHBOURHOOD_DEVIATION) &&
					(std::abs(self->points[i]->y - other->points[j]->y) < MATCH_NEIGHBOURHOOD_DEVIATION)) {
				//				std::cout << "Found match " <<
				//						self->points[i]->x << "," << self->points[i]->y << std::endl;
				self->points[i]->match = other->points[j];
				other->points[i]->match = self->points[j];
			}
		}
	}
}


void DetectLineModuleExt::prepareSegments() {
	if (segmentation != GLUE_POINTS) return;

	Accumulator<DecPoint> & acc = *hough->getAccumulator();

	std::vector<Cell<DecPoint> * >temp;
	temp.resize(acc.size());
	ref(acc.begin(), acc.end(), temp.begin());

	std::cout << "Calculate skew for all points" << std::endl;
	std::for_each(temp.begin(), temp.end(), std::bind1st(std::mem_fun(&DetectLineModuleExt::calculateSkew), this) );

	std::cout << "Find matches for points in different cells" << std::endl;
	combine_pairwise(temp.begin(), temp.end(), mem_fun_bind2 ( &DetectLineModuleExt::findMatches, *this ) );
}

//! Plot the accumulator values as an image
void DetectLineModuleExt::plotAccumulator() {
	Accumulator<DecPoint> & acc = *hough->getAccumulator();
	ASize asize = acc.getSize();
	std::cout << "Create image for accumulator of size " << asize.x << 'x' << asize.y << std::endl;
	CRawImage *a_img = new CRawImage(asize.x, asize.y, 1);
	for (int i = 0; i < asize.x; ++i) {
		for (int j = 0; j < asize.y; ++j) {
			Cell<DecPoint> & c = acc.get(i,j);
			if (c.hits > HIT_THRESHOLD) {
				//					std::cout << "Hit at " << i << ',' << j << " with #hits: " << c.hits << std::endl;
				a_img->setValue(i,j,100);
			}
		}
	}
	a_img->saveBmp("accumulator.bmp");
	delete a_img;
}

/**
 * Consider each cell in the accumulator representing another lines. Print all points that belong to a particular 
 * cell.
 */
void DetectLineModuleExt::writeAssignments(std::string ofilename) {

	std::ofstream ofile;
	ofile.open(ofilename);
	Accumulator<DecPoint> & acc = *hough->getAccumulator();
	ASize asize = acc.getSize();
	std::cout << "Create image for accumulator of size " << asize.x << 'x' << asize.y << std::endl;
	for (int i = 0; i < asize.x; ++i) {
		for (int j = 0; j < asize.y; ++j) {
			Cell<DecPoint> & c = acc.get(i,j);
//			if (c.points.size() > HIT_THRESHOLD) {
				//ofile << c.points.size() << '\t' << (j * asize.x + i) << std::endl;
				for (int p = 0; p < c.points.size(); ++p) {
					ofile << c.points[p]->index << '\t' << c.points[p]->cluster_index << '\t' << 
						c.points[p]->x << '\t' << c.points[p]->y << '\t' << (j * asize.x + i) << std::endl;
				}
//			}
		}
	}
	ofile.close();
}

//! Plot the segments using values in the accumulator
void DetectLineModuleExt::plotSegments() {
	Accumulator<DecPoint> & acc = *hough->getAccumulator();
	ASize asize = acc.getSize();

	int height = image->getheight();
	int width = image->getwidth();
	CRawImage *l_img = new CRawImage(width, height, 1);

	// first print raster, this explains some artifacts which are out there
	for (int i = 0; i < width; ++i) {
		for (int j = 0; j < height; ++j) {
			if ((i % default_patch.width) == 0) l_img->setValue(i,j,20);
			if ((j % default_patch.height) == 0) l_img->setValue(i,j,20);
		}
	}

	std::cout << "Print " << segments.size() << " segments" << std::endl;
	// print segments
	for (int l = 0; l < segments.size(); l++) {
		//std::cout << "Plot line " << segments[l].src->x << "," << segments[l].src->y << " to " <<
		//		segments[l].dest->x << "," << segments[l].dest->y << std::endl;
		l_img->plotLine(segments[l].src->x, segments[l].src->y, segments[l].dest->x, segments[l].dest->y);
	}

	//	std::vector<DecPoint*> & pnts = hough.getPoints();
	//	std::cout << "Number of points: " << pnts.size() << std::endl;
	//	for (int l = 0; l < pnts.size(); l++) {
	//		if (pnts[l]->match != NULL) {
	//			std::cout << "Match at " << pnts[l]->x << "," << pnts[l]->y << std::endl;
	//			l_img->plotCross(pnts[l]->x, pnts[l]->y, 3);
	//		}
	//	}

	l_img->saveBmp("reconstructed.bmp");
	delete l_img;
}


bool DetectLineModuleExt::Stop() {
	if (stop) {
		delete image;
		std::cout << "Gracefully stopped" << std::endl;
		return true;
	}
	return false;
}

