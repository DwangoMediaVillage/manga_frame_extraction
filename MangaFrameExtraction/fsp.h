//
//  fsp.h
//  MangaFrameExtraction
//
//  Created by 山田　祐雅　 on 2015/11/02.
//  Copyright (c) 2015年 山田　祐雅　. All rights reserved.
//
//  reference: http://www.ams.giti.waseda.ac.jp/pdf-files/ishii_dissertation_final.pdf

#ifndef __MangaFrameExtraction__fsp__
#define __MangaFrameExtraction__fsp__

#include <stdio.h>
#include <iostream>
#include <opencv2/core/core.hpp>
#include <opencv2/imgproc/imgproc.hpp>
#include <opencv2/highgui/highgui.hpp>

using namespace std;
using namespace cv;

#define COLOR_BLACK CV_RGB(0, 0, 0)
#define COLOR_WHITE CV_RGB(255, 255, 255)
#define AREA_THRESHOLD 1
#define ADD_PAGEFRAME_WIDTH 20
#define N_BIN 45
#define THETA (180 / N_BIN)
#define BLOCK_SIZE 3
#define CELL_SIZE 1
#define R (CELL_SIZE * (BLOCK_SIZE)*0.5)
#define MARGIN 1
#define NUM_SLC 3
#define NUM_CANDIDATE 10

#define to_dig(rad) ((rad)*180.0 / CV_PI)
#define to_rad(dig) ((dig)*CV_PI / 180.0)

// 画素
typedef struct PixPoint {
    unsigned short x;
    unsigned short y;
#ifdef DEBUG
    operator CvPoint()
    {
        return CvPoint(this->x, this->y);
    }
#endif
} PixPoint;

// 分割線
class SL {
public:
    SL(const SL& sl)
    {
        this->is_horizontal = sl.is_horizontal;
        this->position = sl.position;
        this->theta = sl.theta;
        this->ig = sl.ig;
        this->wpr = sl.wpr;
        this->Hw = sl.Hw;
        this->pixels = sl.pixels;
    }
    SL(bool is_horizontal = true, int position = 0, int theta = 0, double ig = 0.0, double wpr = 0.0, double Hw = 0.0, vector<PixPoint> pixels = vector<PixPoint>(0))
    {
        this->is_horizontal = is_horizontal;
        this->position = position;
        this->theta = theta;
        this->ig = ig;
        this->wpr = wpr;
        this->Hw = Hw;
        this->pixels = pixels;
    }
    ~SL() {}

    SL& operator=(const SL& sl)
    {
        this->is_horizontal = sl.is_horizontal;
        this->position = sl.position;
        this->theta = sl.theta;
        this->ig = sl.ig;
        this->wpr = sl.wpr;
        this->Hw = sl.Hw;
        this->pixels = sl.pixels;
        return (*this);
    }
    bool is_horizontal;
    int position;
    int theta;
    double ig;
    double wpr;
    double Hw;
    vector<PixPoint> pixels;
};

enum Response {
    OK = 0,
    DROP_SLICE_SRC,
    DROP_REMAINED_SRC,
    INVALID_DIGREES,
    INVASION_FRAMES
};

static int separate_count = 0;

class FrameSeparation {
public:
    FrameSeparation(IplImage src /*, IplImage *remained_src*/, string filename, string output_dir, int original_size);
    ~FrameSeparation();

    void save_image(IplImage* img);

    // Center Weighted concentration Gradient Value
    void cwgv();
    // Ditection of Separate Line Candidate for horizontal and vertical direction
    // 分割線候補検出
    void dslc_hv();
    // Ditection of Separate Line Candidate for oblique direction
    // 分割線候補検出
    void dslc_o();

    // コマ内外判定
    int invasion_test(bool is_horizontal, int position, int length, int theta);

    // 傾きのある直線上の画素を走査
    inline void detect_pixels(bool is_horizontal, int position, int length, int theta, vector<PixPoint>* pixels);

    // Separate Line Adaptation Test
    // 分割線適合検査
    void slat();
    bool sl_exists();
    int separation();
    bool is_blank(IplImage* image);

    // 共用の輝度勾配値を計算する
    void calculate_ig();

private:
    void calculate_slc(bool is_horizontal);
    void calculate_wpr(bool is_horizontal);
    void calculate_oblique_slc(bool is_horizontal);
    void calculate_oblique_wpr(bool is_horizontal);

    inline CvPoint calculate_ending_point(bool is_horizontal, int position, int length, int theta)
    {
        return is_horizontal ? CvPoint(length, position + length * cos(to_rad(theta))) : CvPoint(position + length * cos(to_rad(theta)), length);
    }

    // 元画像
    IplImage* src;
    // 分割線で切った画像
    IplImage* slice_src;
    // 分割線で切った残りの画像
    IplImage* remained_src;
    // 作業用画像
    IplImage* proc_img;
    // 二値化画像
    IplImage* bin_img;
    // detect_pixels用
    Mat dp_img;

    // HOG
    vector<Mat> integrals;

    int original_size;

    // for cwgv
    int xi;

    // 分割候補線
    vector<SL> slc[2];
    // x,y軸の各分割線候補の位置
    int sl_position[2];
    // x,y軸の各分割線候補の評価値
    double sl_hw[2];

    SL slice_line;

    FrameSeparation* fs1_recursive;
    FrameSeparation* fs2_recursive;

    string filename;
    string output_dir;

    // 共用の輝度勾配値
    Mat ig_mat;
};

#endif /* defined(__MangaFrameExtraction__fsp__) */
