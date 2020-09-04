#include "CRawImage.h"
#include <cassert>
#include <iostream>

#define RGB_HEADER_SIZE		54
#define PALETTE_SIZE		(256*4)
#define GRAY_HEADER_SIZE	(RGB_HEADER_SIZE+PALETTE_SIZE)
#define BITMAP_SIZE			0

int CRawImage::numSaved = 0;

/**
 * Add a simple ASSERT macro, to see directly where the assert exactly fails.
 */
#define ASSERT(condition) { \
	if(!(condition)){ \
		std::cerr << "ASSERT FAILED: " << #condition << " @ " << __FILE__ << " (" << __LINE__ << ")" << std::endl; \
		assert(condition); \
	} \
	}

/**
 * Use just a raw header so we don't get packed structs etc. By default we will assume an RGB format. Comments are
 * added to the values, because that one makes it very obfuscated where the fields are actually for.
 */
static VALUE_TYPE header[] = {66,77, // magic word for BMP
		RGB_HEADER_SIZE,0,0,0, // size of the BMP file (no data is just the header)
		0,0, // application specific
		0,0, // application specific
		54,0,0,0, // offset to the bitmap
		40,0,0,0, // number of bytes in the DIB header
		0,0,0,0, // width in pixels (18,19,20,21)
		0,0,0,0, // height in pixels
		1,0, // number of color planes, 1 is only legal value
		24,0, // 24-bits (number of bits per pixel), by default RGB with each 8 bits
		0,0,0,0, // BI_RGB, no compression
		0,0,0,0, // size of the raw data in the pixel array including padding, including bitmap
		19,11,0,0, //horizontal resolution
		19,11,0,0, //vertical resolution
		0,0,0,0, //number of colors
		0,0,0,0, // 0 means all colors are important
		0,0,255, //red (start of the bitmap)
		255,255,255, //white pixel
		0,0, // padding for 4-byte alignment
		255,0,0, //blue
		0,255,0, //green
		0,0}; // padding

/**
 * Size will be set automatically and concerns not the number of pixels, but the memory
 * space required, so: width*height*bpp.
 */
CRawImage::CRawImage(int wi, int he, int bpp): width(wi), height(he), bpp(bpp)
{
	size = bpp*width*height;
	data = (VALUE_TYPE*)calloc(size,sizeof(VALUE_TYPE));
	updateHeader();
}

CRawImage::CRawImage(const char *name) {
	data = (VALUE_TYPE*)calloc(PALETTE_SIZE,sizeof(VALUE_TYPE));
	loadBmp(name);
}

CRawImage::CRawImage(const CRawImage & other): width(other.width), height(other.height),
	  size(other.size), bpp(other.bpp) {
	  data = (VALUE_TYPE*)calloc(size,sizeof(VALUE_TYPE));
	  memcpy (data, other.data, other.size);
	  updateHeader();
}

VALUE_TYPE CRawImage::getValue(int x, int y) {
	assert (bpp == 1);
	assert (x >= 0 && x < width);
	assert (y >= 0 && y < height);
	int index = (y*width+x);
	return data[index];
}

Pixel CRawImage::getPixel(int x, int y) {
	assert (bpp == 3);
	assert (x >= 0 && x < width);
	assert (y >= 0 && y < height);
	int base = (y*width+x)*3;
	Pixel pixel;
	pixel.b = data[base+0];
	pixel.r = data[base+1];
	pixel.g = data[base+2];
	return pixel;
}

Pixel CRawImage::getPixel(int x, int y, Patch & patch) {
	assert (x >= 0 && x < patch.width);
	assert (y >= 0 && y < patch.height);
	int base = (y*patch.width+x)*3;
	Pixel pixel;
	pixel.b = patch.data[base+0];
	pixel.r = patch.data[base+1];
	pixel.g = patch.data[base+2];
	return pixel;
}

void CRawImage::setPixel(int x, int y, Patch & patch, Pixel pixel) {
	ASSERT(bpp == 3);
	ASSERT(x >= 0 && x < patch.width);
	ASSERT(y >= 0 && y < patch.height);
	int base = (y*patch.width+x)*bpp;
	patch.data[base+0] = pixel.b;
	patch.data[base+1] = pixel.r;
	patch.data[base+2] = pixel.g;
}

void CRawImage::setValue(int x, int y, Patch & patch, VALUE_TYPE value) {
	ASSERT(bpp == 1);
	ASSERT(x >= 0 && x < patch.width);
	ASSERT(y >= 0 && y < patch.height);
	int base = (y*patch.width+x)*bpp;
	patch.data[base] = value;
}

void CRawImage::setValue(int x, int y, VALUE_TYPE value) {
	assert (bpp == 1);
	assert (x >= 0 && x < width);
	assert (y >= 0 && y < height);
	int index = (y*width+x);
	data[index] = value;
}

/**
 * Get patch of size "patch.width * patch.length". This is more "bewerkelijk" than seems at first sight because the
 * array is stored in a linear fashion, to get a "square" of data out of it, you will need to jump all over the place.
 * That's why we now just check every individual element.
 */
void CRawImage::getPatch(int x, int y, Patch & patch) {
//	std::cout << "Get patch at " << x << ',' << y << std::endl;
	assert (x >= 0 && x < width);
	assert (y >= 0 && y < height);
	patch.data = (VALUE_TYPE*)calloc(patch.width*patch.height*bpp, sizeof(VALUE_TYPE));

	for (int i = 0; i < patch.width; ++i) {
		for (int j = 0; j < patch.height; ++j) {
			int base = ((y+j)*width+x+i)*bpp;
			if (bpp == 3) {
				Pixel pixel;
				pixel.b = data[base+0];
				pixel.r = data[base+1];
				pixel.g = data[base+2];
				setPixel(i,j, patch, pixel);
			} else {
				setValue(i,j,patch,data[base]);
			}
		}
	}
}

CRawImage* CRawImage::patch2Img(Patch &patch) {
	CRawImage *img = new CRawImage(patch.width, patch.height, 3);
	memcpy(img->data, patch.data, patch.width*patch.height*3);
	return img;
}


void CRawImage::refresh() {
	size = bpp*width*height;
	if (data != NULL) free(data);
	data = (VALUE_TYPE*)calloc(size,sizeof(VALUE_TYPE));
	std::cout << "Set size to " << width << '*' << height << '*' << bpp << std::endl;
	updateHeader();
}

void CRawImage::readHeader(VALUE_TYPE *hdr) {
	std::cout << "Magic header is " << (int)hdr[0] << ' ' << (int)hdr[1] << std::endl;
	std::cout << "Header size " << (int)hdr[2] << std::endl;
	std::cout << "Offset to bitmap " << (int)hdr[10] << std::endl;
	std::cout << "DIB " << (int)hdr[14] << std::endl;
	width = hdr[18] + (hdr[19] << 8);
	height = hdr[22] + (hdr[23] << 8);
	int bits_per_pixel = hdr[28];
	bpp = bits_per_pixel >> 3;
	std::cout << "Dimensions: " << width << '*' << height << " and bits per pixel: " << bits_per_pixel << std::endl;
}

void CRawImage::updateHeader() {
	assert (width > 0);
	int s;
	if (bpp == 3)
		s = height*width*bpp+RGB_HEADER_SIZE;
	else
		s = height*width*bpp+GRAY_HEADER_SIZE;

	header[2] = (VALUE_TYPE)s;
	header[3] = (VALUE_TYPE)(s >> 8);
	header[4] = (VALUE_TYPE)(s >> 16);

	if (bpp == 3) {
		header[10] = (VALUE_TYPE)RGB_HEADER_SIZE;
		header[11] = (VALUE_TYPE)(RGB_HEADER_SIZE >> 8);
	} else {
		// we make it plus 4*256 values for color pallete in grayscale case
		header[10] = (VALUE_TYPE)GRAY_HEADER_SIZE;
		header[11] = (VALUE_TYPE)(GRAY_HEADER_SIZE >> 8);
	}

	header[18] = width%256;
	header[19] = width/256;
//		header[18] = 255;
//		header[19] = 255;
	// image height 22,23, 24, 25
	header[22] = height%256;
	header[23] = height/256;
	// planes, 26,27
	// bits (not bytes) per pixel 28,29
	header[28] = (bpp*8)%256;
	header[29] = (bpp*8)/256;
	// compression 30,31,32,33
	//  size of the raw data, 34,35,36,37
	s = height*width*bpp;
	header[34] = s%256;
	header[35] = s >> 8;
	header[36] = s >> 16;
	// horizontal resolution
	if (bpp == 1) {
//		header[38] = 18;
//		header[42] = 18;
		for (int i = 44; i < RGB_HEADER_SIZE; ++i) {
			header[i] = 0;
		}
		header[47] = 1;
		header[51] = 1;
	}
}

int CRawImage::getSaveNumber()
{
	char name[100];
	FILE* file = NULL;
	do{
		sprintf(name,"%04i.bmp",numSaved++);
		file = fopen(name,"r");
	}
	while (file != NULL);
	numSaved--;
	return numSaved;
}

CRawImage::~CRawImage()
{
	if (data != NULL) {
		free(data);
	}
}

void CRawImage::makeMonochrome() {
	assert (bpp == 3);
	assert (data != NULL);
	assert (width > 0);
	printf("Size is %i\n", size);
	swap();

	VALUE_TYPE* newData = (VALUE_TYPE*)calloc(width*height,sizeof(VALUE_TYPE));
	for (int i = 0; i < width*height; ++i) {
		int temp = data[i*3]*30 + data[i*3+1]*59 + data[i*3+2]*11;
		newData[i] = temp / (100);
	}
	setbpp(1);
	memcpy(data,newData,size);
	if (newData != NULL) free(newData);
}

void CRawImage::makeMonochrome(CRawImage *result) {
	if (result == NULL) {
		fprintf(stderr, "Resulting image should already be constructed\n");
		return;
	}
	if (result->bpp == 3) result->makeMonochrome();

	if (bpp == 3) {
		for (int i = 0; i < width*height; ++i) {
			int temp = data[i*3]*30 + data[i*3+1]*59 + data[i*3+2]*11;
			result->data[i] = temp / (100);
		}
	} else {
		fprintf(stderr, "Source image is already monochrome, you could've used the copy constructor\n");
		for (int i = 0; i < width*height; ++i) {
			result->data[i] = data[i];
		}
	}
}

/**
 * Swap the first and third bytes of an YUV or RGB array
 */
void CRawImage::swap()
{
	if (bpp != 3) {
		return;
	}
	if (data == NULL) return;
	VALUE_TYPE* newData = (unsigned char*)calloc(size,sizeof(VALUE_TYPE));
	int span = width*bpp;
	for (int j = 0;j<height;j++){
		memcpy(&newData[span*j],&data[span*(height-1-j)],span);
		for (int i = 0;i<width;i++){
			VALUE_TYPE a = newData[(width*j+i)*3];
			newData[(width*j+i)*3] = newData[(width*j+i)*3+2];
			newData[(width*j+i)*3+2] = a;
		}
	}
	memcpy(data,newData,size);
	if (newData != NULL) free(newData);
}

void CRawImage::saveBmp(const char* inName)
{
	updateHeader();
//	std::cout << __func__ << ": save" << std::endl;
	FILE* file = fopen(inName, "wb");
	if (file == NULL) {
		fprintf(stderr, "Could not open file %s", inName);
		return;
	}
//	std::cout << __func__ << ": save2" << std::endl;
	swap();
//	if (bpp == 1)
//		fwrite(bw_header,RGB_HEADER_SIZE,1,file);
//	else
	fwrite(header,RGB_HEADER_SIZE,1,file);
//#ifdef READ_WRITE_PALLETTE
	if (bpp == 1) {
		VALUE_TYPE *temp_palette = (VALUE_TYPE*)calloc(PALETTE_SIZE, sizeof(VALUE_TYPE));
		// you'll need a color palette for grayscale images.
		for (int i = 0; i < PALETTE_SIZE / 4; ++i) {
			temp_palette[i*4] = temp_palette[i*4+1] = temp_palette[i*4+2] = i;
			temp_palette[i*4+3] = 0;
		}
		fwrite(temp_palette, PALETTE_SIZE, 1, file);
		if (temp_palette != NULL) free(temp_palette);
	}
//#endif
	fwrite(data,size,1,file);
	swap();
	fclose(file);
	std::cout << __func__ << ": saved \"" << inName << "\"" << std::endl;
}

void CRawImage::saveNumberedBmp(const char *name, bool increment)
{
	if (increment) numSaved++;
	char nameext[100];
	sprintf(nameext,"%s%04i.bmp",name,numSaved);
	saveBmp(nameext);
}

/**
 * Todo: Actually, rewrite this entire class. Make sure images created by convert etc. can be created correctly.
 * Observe that the data array is first used for the header, and then overwritten by the data. So, header information
 * is only stored in the form of fields in this class.
 */
bool CRawImage::loadBmp(const char* inName)
{
	printf("Open file %s\n", inName);

	// we assume that the image to be loaded has 3 bytes per pixel, if
	// currently at bpp = 1, we adjust
//	if (bpp == 1) setbpp(3);

	FILE* file = fopen(inName,"rb");
	if (file!=NULL)
	{
		size_t n = fread(data,RGB_HEADER_SIZE,1,file);
		if (n == 0) return false;
		readHeader(data);
		// can safely be called, because header is not preserved afterwards
		refresh();

//#ifdef READ_WRITE_PALLETTE
		if (bpp == 1) { // there is a "color" palette if the image is grayscale
			n = fread(data,PALETTE_SIZE,1,file);
			if (n == 0) return false;
		}
//#endif
		n = fread(data,size,1,file);
		if (n == 0) return false;
		fclose(file);
//		swap();
		return true;
	}
	return false;
}

void CRawImage::plotCenter()
{
	int centerWidth = 20;
	VALUE_TYPE color[] = {255,150,150};
	for (int i = -centerWidth;i<centerWidth;i++){
		for (int j =0;j<3;j++){
			data[(width*(height/2+i)+width/2-centerWidth)*3+j] = color[j];
			data[(width*(height/2+i)+width/2+centerWidth)*3+j] = color[j];
			data[(width*(height/2-centerWidth)+width/2+i)*3+j] = color[j];
			data[(width*(height/2+centerWidth)+width/2+i)*3+j] = color[j];
		}
	}
}

void CRawImage::plotLine(int x,int y) {
	int base;
	if (y < 0 || y > height-1) y = height/2;
	if (x < 0 || x > width-1) x = width/2;
	for(int i=0; i < width;i++) {
		if (i == width/2) i++;
		base = (width*y+i)*3;
		data[base+0] = 255;
		data[base+1] = 0;
		data[base+2] = 255;
	}

	for(int j=0;j<height;j++) {
		const int bidx = ((width*j)+x)*3;
		if (j == height/2) j++;
		data[bidx+0] = 255;
		data[bidx+1] = 255;
		data[bidx+2] = 0;
	}
}

void CRawImage::plotLine(int x0, int y0, int x1, int y1) {
	assert(bpp == 1); // for now
	assert (x0 >= 0 && x0 < width);
	assert (x1 >= 0 && x1 < width);
	assert (y0 >= 0 && y0 < height);
	assert (y1 >= 0 && y1 < height);
	int xr = abs(x0-x1);
	int yr = abs(y0-y1);
//	std::cout << "Plot line from [" << x0 << "," << y0 << "] to [" << x1 << "," << y1 << "]" << std::endl;
	if (xr > yr) {
//		std::cout << "Adjust x with 1, and adjust y correspondingly slower" << std::endl;
		int xstep = x0 < x1 ? 1 : -1;
		int ystep = y0 < y1 ? 1 : -1;
		for (int i = x0; i != x1; i+=xstep) { // loop over slowest moving index
			int j = y0 + xstep*ystep*((i-x0) * yr) / xr;
//			std::cout << " [" << i << "," << j << "] ";
			setValue(i,j,255);
		}
	} else {
//		std::cout << "Adjust y with 1, and adjust x correspondingly slower" << std::endl;
		int xstep = x0 < x1 ? 1 : -1;
		int ystep = y0 < y1 ? 1 : -1;
//		std::cout << "y will be " << ((ystep < 0) ? "decremented" : "incremented") << std::endl;
//		std::cout << "x will be " << ((xstep < 0) ? "decremented" : "incremented") << std::endl;
		for (int j = y0; j != y1; j+=ystep) { // loop over slowest moving index
			int i = x0 + xstep*ystep*((j-y0) * xr) / yr;
//			std::cout << " [" << i << "," << j << "] ";
			setValue(i,j,255);
		}
	}
//	std::cout << std::endl;
}

/**
 * Plotting a (little) cross at position (i,j). The size of the cross is in pixels and is "just one leg" of it. The
 * entire cross with size=3, becomes 7 pixels wide. The number of bytes per pixel, bpp can be 3 or 1.
 */
void CRawImage::plotCross(int i,int j,int size) {
	int cross = 4;
	for (int di = -size; di < size; ++di) {
		int dii = i+di;
		if (dii < 0) continue;
		if (dii >= width) continue;
		if (bpp == 3) {
			data[dii*bpp+j*width*bpp+0] = 0;
			data[dii*bpp+j*width*bpp+1] = 0;
			data[dii*bpp+j*width*bpp+2] = 255;
		} else {
			data[dii*bpp+j*width*bpp+0] = 255;
		}
	}
	for (int dj = -size; dj < size; ++dj) {
		int djj = j+dj;
		if (djj < 0) continue;
		if (djj >= height) continue;
		if (bpp == 3) {
			data[i*bpp+djj*width*bpp+0] = 0;
			data[i*bpp+djj*width*bpp+1] = 0;
			data[i*bpp+djj*width*bpp+2] = 255;
		} else {
			data[i*bpp+djj*width*bpp+0] = 255;
		}
	}
}

/** pocita jas obrazku:
 *  upperHalf == true, pocita se jen z horni poloviny obrazku
 *  upperHalf == false, pocita jen ze spodni poloviny obrazku
 */
double CRawImage::getOverallBrightness(bool upperHalf) {
	int step = 5;
	int sum,num,satMax,satMin,pos;
	sum=num=satMax=satMin=0;
	int limit = 0;
	if (upperHalf) limit = 0; else limit=height/2;
	for (int i = limit;i<height/2+limit;i+=step){
		for (int j = 0;j<width;j+=step){
			pos = (i*width+j)*bpp;
			if (data[pos] >= 250 && data[pos+1] >=250 && data[pos+2] >= 250) satMax++;
			if (data[pos] <= 25 && data[pos+1] <=25 && data[pos+2] <= 25) satMin++;
			sum+=data[pos] + data[pos+1] + data[pos+2];
			num++;
		}
	}
	return (sum/num/bpp) + (satMax-satMin)*100.0/num;
}



