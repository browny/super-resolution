
#include <iostream>
#include <string>
#include <math.h>
using namespace std;

// CImag
#include "CImg.h"
using cimg_library::CImg;
using namespace cimg_library;

// Matrix Library
#include "2dmem_template_class.h"

// Construct Weight Table
float WeightTable[256];

void nlmZoom2(CImg<unsigned char>*, CImg<unsigned char>*);
void zoom(CImg<unsigned char>*, CImg<unsigned char>*, int scale);

int main(int argc, char* argv[])
{
	int index;
	for (int i = 0; i < 256; i++) {
		index = i;
		WeightTable[i] = exp(-index * 0.2);
	}

	CImg<unsigned char> inImg;
	CImg<unsigned char> outImg;
	CImg<unsigned char> bicubicImg;

	const char* src = "raw.bmp";
	inImg.load(src);
	bicubicImg.load(src);

	zoom(&inImg, &outImg, 4); // if you want x4 => zoom(&inImg, &outImg, 4);

	const char* result = "out.bmp";
	outImg.save(result);

	bicubicImg.resize(bicubicImg.width*4, bicubicImg.height*4, 1, 3, 5, -1, false);
	bicubicImg.save("bicubic.bmp");

	return 0;
}

void nlmZoom2(CImg<unsigned char>* pInImg, CImg<unsigned char>* pOutImg) {

	const int srch_range = 1;
	CImg<unsigned char> InImg;
	CImg<unsigned char> BicubicImg;
	CImg<unsigned char> OutImg;

	InImg = *pInImg;

	// bicubic x2
	InImg.resize(InImg.width*2, InImg.height*2, 1, 3, 5, -1, false);

	// bicubic x4
	BicubicImg = InImg;
	BicubicImg.resize(BicubicImg.width*2, BicubicImg.height*2, 1, 3, 3, -1, false);

	const int outSize = (BicubicImg.height)*(BicubicImg.width);

	float* Intensity2D = new float [outSize];
	float* Weight2D = new float [outSize];

	for (int ch = 0; ch < 3; ch++) {

		CImg<unsigned char> InImgChannel;		
		InImgChannel = InImg.get_channel(ch);
		CImg<unsigned char> BiImgChannel;		
		BiImgChannel = BicubicImg.get_channel(ch);

		memset(Intensity2D, 0, sizeof(float)*outSize);
		memset(Weight2D, 0, sizeof(float)*outSize);

		int encls_candi;
		float weight;
		int refPoint;

		float Nei_TL_row, NeiLR_start_row;
		float Nei_TL_col, NeiLR_start_col;
		float Nei_BR_row, NeiLR_end_row  ;
		float Nei_BR_col, NeiLR_end_col  ;
		
		for (int i = 2; i <= InImg.height * 2 - 2; i += 2) {
			for (int j = 2; j <= InImg.width * 2 - 2; j += 2) {
				
				// On
				refPoint = BiImgChannel(j, i, 0, 0);
				
				Nei_TL_row = i - srch_range;
				Nei_TL_col = j - srch_range;
				Nei_BR_row = i + srch_range;
				Nei_BR_col = j + srch_range;

				NeiLR_start_row = Nei_TL_row/2;
				NeiLR_start_col = Nei_TL_col/2;
				NeiLR_end_row   = Nei_BR_row/2;
				NeiLR_end_col   = Nei_BR_col/2;
				
				for (int k = NeiLR_start_row; k <= NeiLR_end_row; k++) {
					for (int l = NeiLR_start_col; l <= NeiLR_end_col; l++) {

						encls_candi = InImgChannel(l, k, 0, 0);
						weight = WeightTable[abs(refPoint - encls_candi)];

						for (int paste_i = i - 1; paste_i <= i + 1; paste_i++) {
							for (int paste_j = j - 1; paste_j <= j + 1; paste_j++) {
								
								Intensity2D[(paste_i - 1) * (InImg.width * 2)
										+ (paste_j - 1)] += weight
										* encls_candi;

								Weight2D[(paste_i - 1) * (InImg.width * 2)
										+ (paste_j - 1)] += weight;

							}
						}
					}
				}
			
				// Horizontal
				refPoint = BiImgChannel(j, i + 1, 0, 0);

				Nei_TL_row = i+1 - srch_range;
				Nei_TL_col = j - srch_range;
				Nei_BR_row = i+1 + srch_range;
				Nei_BR_col = j + srch_range;

				NeiLR_start_row = 1 + ceil((Nei_TL_row-1)/2);
				NeiLR_start_col = 1 + ceil((Nei_TL_col-1)/2);
				NeiLR_end_row   = 1 + floor((Nei_BR_row-1)/2);
				NeiLR_end_col   = 1 + floor((Nei_BR_col-1)/2);

				for (int k = NeiLR_start_row; k <= NeiLR_end_row; k++) {
					for (int l = NeiLR_start_col; l <= NeiLR_end_col; l++) {
						encls_candi = InImgChannel(l, k, 0, 0);
						weight = WeightTable[abs(refPoint - encls_candi)];

						for (int paste_i = i - 1; paste_i <= i + 1; paste_i++) {
							for (int paste_j = (j + 1) - 1; paste_j <= (j + 1)
									+ 1; paste_j++) {
								Intensity2D[(paste_i - 1) * (InImg.width * 2)
										+ (paste_j - 1)] += weight
										* encls_candi;
								Weight2D[(paste_i - 1) * (InImg.width * 2)
										+ (paste_j - 1)] += weight;

							}
						}

					}
				}

				// Vertical
				refPoint = BiImgChannel(j+1, i, 0, 0);

				Nei_TL_row = i - srch_range;
				Nei_TL_col = j+1 - srch_range;
				Nei_BR_row = i + srch_range;
				Nei_BR_col = j+1 + srch_range;

				NeiLR_start_row = 1 + ceil((Nei_TL_row-1)/2);
				NeiLR_start_col = 1 + ceil((Nei_TL_col-1)/2);
				NeiLR_end_row   = 1 + floor((Nei_BR_row-1)/2);
				NeiLR_end_col   = 1 + floor((Nei_BR_col-1)/2);


				for (int k = NeiLR_start_row; k <= NeiLR_end_row; k++) {
					for (int l = NeiLR_start_col; l <= NeiLR_end_col; l++) {
						encls_candi = InImgChannel(l, k, 0, 0);
						weight = WeightTable[abs(refPoint - encls_candi)];

						for (int paste_i = (i + 1) - 1; paste_i <= (i + 1) + 1; paste_i++) {
							for (int paste_j = j - 1; paste_j <= j + 1; paste_j++) {
								Intensity2D[(paste_i - 1) * (InImg.width * 2)
										+ (paste_j - 1)] += weight
										* encls_candi;
								Weight2D[(paste_i - 1) * (InImg.width * 2)
										+ (paste_j - 1)] += weight;

							}
						}
					}
				}

				// Cross
				
				refPoint = BiImgChannel(j+1, i+1, 0, 0);

				Nei_TL_row = i+1 - srch_range;
				Nei_TL_col = j+1 - srch_range;
				Nei_BR_row = i+1 + srch_range;
				Nei_BR_col = j+1 + srch_range;

				NeiLR_start_row = Nei_TL_row/2;
				NeiLR_start_col = Nei_TL_col/2;
				NeiLR_end_row   = Nei_BR_row/2;
				NeiLR_end_col   = Nei_BR_col/2;
				

				for (int k = NeiLR_start_row; k <= NeiLR_end_row; k++) {
					for (int l = NeiLR_start_col; l <= NeiLR_end_col; l++) {
						encls_candi = InImgChannel(l, k, 0, 0);
						weight = WeightTable[abs(refPoint - encls_candi)];

						for (int paste_i = (i + 1) - 1; paste_i <= (i + 1) + 1; paste_i++) {
							for (int paste_j = (j + 1) - 1; paste_j <= (j + 1)
									+ 1; paste_j++) {
								Intensity2D[(paste_i - 1) * (InImg.width * 2)
										+ (paste_j - 1)] += weight
										* encls_candi;
								Weight2D[(paste_i - 1) * (InImg.width * 2)
										+ (paste_j - 1)] += weight;

							}
						}
					}
				}
			}
		}		

		for (int i = 0; i < InImg.height * 2; i++)
			for (int j = 0; j < InImg.width * 2; j++)
				BicubicImg(j, i, 0, ch) = (unsigned char)
				((Intensity2D[i * InImg.width * 2 + j] / Weight2D[i * InImg.width * 2 + j]));

	}

	delete[] Intensity2D;
	delete[] Weight2D;

	BicubicImg.resize(BicubicImg.width/2, BicubicImg.height/2, 1, 3, 5, -1, false);
	*pOutImg = BicubicImg;

	
}

void zoom(CImg<unsigned char>* InImg, CImg<unsigned char>* OutImg, int scale) {

	const int numLoop = scale / 2;

	/*long cTick, nTick;
	 cTick = GetTickCount();*/

	for (int l = 0; l < numLoop; l++) {

		// Reset
		if (l == 0) {

			OutImg->assign(InImg->width * 2, InImg->height * 2, 1, 3);
			nlmZoom2(InImg, OutImg);

		} else {

			InImg->clear();
			*InImg = *OutImg;
			OutImg->clear();

			OutImg->assign(OutImg->width * 2, OutImg->height * 2, 1, 3);
			nlmZoom2(InImg, OutImg);
		}

		/*nTick = GetTickCount();
		 cout << "executed in" << nTick - cTick << "ms \n";*/
	}
}
