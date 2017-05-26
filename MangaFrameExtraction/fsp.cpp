//
//  fsp.cpp
//  MangaFrameExtraction
//
//  Created by 山田　祐雅　 on 2015/11/02.
//  Copyright (c) 2015年 山田　祐雅　. All rights reserved.
//

#include "fsp.h"

FrameSeparation::FrameSeparation(IplImage src, string filename, string output_dir, int original_size, PixPoint rel_original_point)
{

    this->src = cvCloneImage(&src);
    //    this->remained_src = cvCloneImage(remained_src);
    this->proc_img = cvCloneImage(&src);
    this->bin_img = cvCloneImage(&src);
    //    this->separate_count = separate_count;
    this->original_size = original_size;
    this->rel_original_point = rel_original_point;
    this->filename = filename;
    this->output_dir = output_dir;
    this->dp_img = cvarrToMat(this->src);
    // 入力画像が真っ白or非常に小さく、分割する必要があるかどうかを確認
    if (!is_blank(this->src)) {
        calculate_ig();

        // 最大長の2%分を走査から外す
        xi = MIN(src.width, src.height) * 0.02;

#ifdef DEBUG
        cout << "===" << separate_count << "===" << endl;
#endif
        cwgv();
        dslc_hv();

        slat();
        if (!sl_exists()) {
            // 斜めのコマを考慮する
            dslc_o();
            slat();
            if (!sl_exists()) {
                if (!is_blank(&src)) {
                    save_image(&src);
                }
            }
            else {
                switch (separation()) {
                case DROP_SLICE_SRC:
                    if (this->remained_src->width * this->remained_src->height >= this->src->width * this->src->height * 0.95) {
                        save_image(this->src);
                        break;
                    }
                    fs1_recursive = new FrameSeparation(*this->remained_src, filename, output_dir, original_size, rel_remained_point);
                    delete fs1_recursive;
                    break;

                case DROP_REMAINED_SRC:
                    if (this->slice_src->width * this->slice_src->height >= this->src->width * this->src->height * 0.95) {
                        save_image(this->src);
                        break;
                    }
                    fs1_recursive = new FrameSeparation(*this->slice_src, filename, output_dir, original_size, rel_slice_point);
                    delete fs1_recursive;
                    break;
                case OK:
                    fs1_recursive = new FrameSeparation(*this->slice_src, filename, output_dir, original_size, rel_slice_point);
                    fs2_recursive = new FrameSeparation(*this->remained_src, filename, output_dir, original_size, rel_remained_point);
                    delete fs1_recursive;
                    delete fs2_recursive;
                    break;
                default:
                    break;
                }
            }
        }
        else {
            switch (separation()) {
            case DROP_SLICE_SRC:
                fs1_recursive = new FrameSeparation(*this->remained_src, filename, output_dir, original_size, rel_remained_point);
                delete fs1_recursive;
                break;

            case DROP_REMAINED_SRC:
                fs1_recursive = new FrameSeparation(*this->slice_src, filename, output_dir, original_size, rel_slice_point);
                delete fs1_recursive;
                break;
            case OK:
                fs1_recursive = new FrameSeparation(*this->slice_src, filename, output_dir, original_size, rel_slice_point);
                fs2_recursive = new FrameSeparation(*this->remained_src, filename, output_dir, original_size, rel_remained_point);
                delete fs1_recursive;
                delete fs2_recursive;
                break;
            default:
                break;
            }
        }
    }
}

FrameSeparation::~FrameSeparation()
{
    cvReleaseImage(&proc_img);
    cvReleaseImage(&src);
    cvReleaseImage(&bin_img);
}

void FrameSeparation::save_image(IplImage* img)
{
    cout << "panel-region" << " ";
    cout << separate_count << " ";
    cout << rel_original_point.x << " ";
    cout << rel_original_point.y << " ";
    cout << img->width << " ";
    cout << img->height << " ";
    cout << endl;

    stringstream fn;
    fn << output_dir << filename << "_" << separate_count << ".jpg";
    cvSaveImage(fn.str().c_str(), img);
    separate_count++;
}

void FrameSeparation::cwgv()
{
    // 2値化
    IplImage* binary = cvCloneImage(src);
    cvThreshold(binary, binary, 120, 255, CV_THRESH_BINARY);
    this->bin_img = cvCloneImage(binary);
    this->proc_img = cvCloneImage(binary);
#ifdef DEBUG
    cvShowImage("[ cwgv ] bin_img", bin_img);
    cv::waitKey(0);
#endif
    cvReleaseImage(&binary);
}

void FrameSeparation::dslc_hv()
{
    // y軸の走査
    calculate_slc(true);
    calculate_wpr(true);
    // x軸の走査
    calculate_slc(false);
    calculate_wpr(false);
}

void FrameSeparation::dslc_o()
{
    // y軸の走査
    calculate_oblique_slc(true);
    //    calculate_oblique_wpr(true);
    // x軸の走査
    calculate_oblique_slc(false);
    //    calculate_oblique_wpr(false);
}

int FrameSeparation::invasion_test(bool is_horizontal, int position, int length, int theta)
{
    vector<PixPoint> pixels;
    try {
        detect_pixels(is_horizontal, position, length, theta, &pixels);
    }
    catch (Response err) {
        if (err == INVALID_DIGREES) {
#ifdef DEBUG
            cout << "invalid digree" << endl;
#endif
            return err;
        }
    }

    bool is_left_black = false, is_right_black = false;
    int width = theta == 90 ? 2 : 3;
    int count = 0;
    int count_l = 0, count_r = 0;
    for (int d = 0; d < pixels.size(); d++) {
        if (pixels[d].x + width >= src->width || pixels[d].x - width <= 0 || pixels[d].y + width >= src->height || pixels[d].y - width <= 0) {
            continue;
        }
        is_left_black = false, is_right_black = false;
        if (is_horizontal) {
            for (int i = 1; i <= width; i++) {
                is_left_black = (uchar) this->bin_img->imageData[(pixels[d].y + i) * this->bin_img->widthStep + pixels[d].x * this->bin_img->nChannels] < 127 ? true : is_left_black;
                is_right_black = (uchar) this->bin_img->imageData[(pixels[d].y - i) * this->bin_img->widthStep + pixels[d].x * this->bin_img->nChannels] < 127 ? true : is_right_black;
            }
        }
        else {
            for (int i = 1; i <= width; i++) {
                is_left_black = (uchar) this->bin_img->imageData[pixels[d].y * this->bin_img->widthStep + (pixels[d].x - i) * this->bin_img->nChannels] < 127 ? true : is_left_black;
                is_right_black = (uchar) this->bin_img->imageData[pixels[d].y * this->bin_img->widthStep + (pixels[d].x + i) * this->bin_img->nChannels] < 127 ? true : is_right_black;
            }
        }
        if (is_left_black && !is_right_black) {
            count_l++;
        }
        if (!is_left_black && is_right_black) {
            count_r++;
        }
    }
    count = MAX(count_l, count_r);
#ifdef DEBUG
    cout << "is_horizontal:" << is_horizontal << ", position:" << position << ", theta:" << theta << ", invasion test:" << count / (double)length << endl;
#endif
    return count / (double)length > (theta == 90 ? 0.55 : 0.4) ? OK : INVASION_FRAMES;
}

inline void FrameSeparation::detect_pixels(bool is_horizontal, int position, int length, int theta, vector<PixPoint>* pixels)
{
    if (theta > 135 || theta < 45) {
        throw(INVALID_DIGREES);
    }
    int x0 = is_horizontal ? 0 : position;
    int y0 = is_horizontal ? position : 0;

    CvPoint point = calculate_ending_point(is_horizontal, position, length, theta);

    if (point.x <= 0 || point.y <= 0) {
        throw(INVALID_DIGREES);
    }

    LineIterator it(dp_img, CvPoint(x0, y0), cvPoint(point.x, point.y), 8);

    pixels->resize(it.count);
    for (int i = 0; i < it.count; i++, ++it) {
        pixels->at(i).x = it.pos().x;
        pixels->at(i).y = it.pos().y;
    }
}

void FrameSeparation::slat()
{
    // 定数
    int delta = 40;
    double rho = 0.20;
    int N = 7;
    int M = 3;
    if (slc[0].size() > src->width) {
        delta = 35;
        rho = 0.30;
        N = 20;
        M = 5;
    }
    if (N <= M + 2) {
        throw "SLAT(): n or m doesnt match the criteria";
    }

    slice_line = SL();

    // 最終的な分割線候補
    vector<SL> slc_final[2];
    int length, ep;
    // 分割線候補から、x,y軸それぞれ最終的な分割線候補を決定する
    for (int i = 0; i < 2; i++) {
#ifdef DEBUG
        cout << (i == 0 ? "x axis" : "y axis") << endl;
#endif
        length = (bool)i ? src->height : src->width;

        // 分割線候補のHwを計算する
        for (int j = xi; j < slc[i].size() - xi; j++) {
            if (slc[i].size() > (i ? src->width : src->height)) {
                slc[i].at(j).Hw = slc[i].at(j).ig;
            }
            else {
                slc[i].at(j).Hw = slc[i].at(j).ig * slc[i].at(j).wpr;
            }
        }
        // 上位NUM_CANDIDATE件を最終的な分割線候補とする
        for (int count = 0; count < NUM_CANDIDATE; count++) {
            double max = 0.0;
            slc_final[i].push_back(SL());
            for (int j = xi; j < slc[i].size() - xi; j++) {
                if (slc[i].at(j).pixels.size() == 0 || slc[i].at(j).position <= xi) {
                    continue;
                }
                // 四隅xiピクセル(から始まる|で終わる）分割候補線は走査から除外
                ep = (bool)i ? slc[i].at(j).pixels.at(slc[i].at(j).pixels.size() - 1).y : slc[i].at(j).pixels.at(slc[i].at(j).pixels.size() - 1).x;
                if (abs(length - slc[i].at(j).position) <= xi || abs(length - ep) <= xi) {
                    continue;
                }
                if (count == 0) {
                    if (slc_final[i].at(count).Hw < slc[i].at(j).Hw) {
                        slc_final[i].data()[count] = slc[i].at(j);
                    }
                }
                else {
                    if (max < slc[i].at(j).Hw) {
                        if (slc_final[i].at(count - 1).Hw > slc[i].at(j).Hw) {
                            max = slc[i].at(j).Hw;
                            slc_final[i].data()[count] = slc[i].at(j);
                        }
                    }
                }
            }
            // 見つからなかった場合
            if (slc_final[i].at(count).position == 0) {
                continue;
            }
#ifdef DEBUG
            cout << "(" << slc_final[i].at(count).position << ", " << slc_final[i].at(count).Hw << ")" << endl;
            if (i == 0) {
                cvLine(proc_img, cvPoint(slc_final[i].at(count).position, 0), (CvPoint)slc_final[i].at(count).pixels.at(slc_final[i].at(count).pixels.size() - 1), CV_RGB(255, 0, 0));
            }
            else {
                cvLine(proc_img, cvPoint(0, slc_final[i].at(count).position), (CvPoint)slc_final[i].at(count).pixels.at(slc_final[i].at(count).pixels.size() - 1), CV_RGB(255, 0, 0));
            }
#endif
        }
    }
#ifdef DEBUG
    cvShowImage("[ slat ] slice line", proc_img);
    cv::waitKey(0);
#endif

    // x,y軸のどちらの分割線が、正しい・分割すべき線かを決定する
    // 評価方法は、分割線上輝度勾配分布による
    int igtest_pass_count_tmp;
    double total_test_retio[2][NUM_CANDIDATE];
    int total_test_pass_count[2][NUM_CANDIDATE];
    double axis_test_retio[2];
    int axis_test_retio_index[2], direction = 0;
    vector<PixPoint>* pixels;
    for (int axis = 0; axis < 2; axis++) {
        for (int candidate = 0; candidate < NUM_CANDIDATE; candidate++) {
            pixels = &slc_final[axis].at(candidate).pixels;
            direction = (int)pixels->size();
            int* igh = new int[direction];
            for (int d = 0; d < direction; d++) {
                if ((bool)axis) {
                    // 横に分割
                    if (ig_mat.at<Vec3b>(pixels->at(d).y, pixels->at(d).x)[1] >= 90) {
                        igh[d] = 180 - ig_mat.at<Vec3b>(pixels->at(d).y, pixels->at(d).x)[1];
                    }
                    else {
                        igh[d] = ig_mat.at<Vec3b>(pixels->at(d).y, pixels->at(d).x)[1];
                    }
                }
                else {
                    // 縦に分割
                    if (ig_mat.at<Vec3b>(pixels->at(d).y, pixels->at(d).x)[1] >= 90) {
                        igh[d] = -90 + ig_mat.at<Vec3b>(pixels->at(d).y, pixels->at(d).x)[1];
                    }
                    else {
                        igh[d] = 90 - ig_mat.at<Vec3b>(pixels->at(d).y, pixels->at(d).x)[1];
                    }
                }
            }
            igtest_pass_count_tmp = 0;
            // 90度の場合、最初と最後の余白（こう配が0）は計算から除外する
            int left_marin_count = 0, right_margin_count = direction - 1;
            if (slc_final[axis].at(candidate).theta == 90) {
                for (int i = (BLOCK_SIZE - 1) / 2; i < direction - (BLOCK_SIZE - 1) / 2;
                     i++) {
                    if (igh[i] == 0) {
                        left_marin_count = i;
                    }
                    else {
                        break;
                    }
                }
                for (int i = direction - (BLOCK_SIZE - 1) / 2 - 1;
                     i >= (BLOCK_SIZE - 1) / 2; i--) {
                    if (igh[i] == 0) {
                        right_margin_count = i;
                    }
                    else {
                        break;
                    }
                }
            }

            int new_count = right_margin_count - left_marin_count;
            for (int n = 0; n < N; n++) {
                int count = 0;
                for (int i = n * (new_count / N) + left_marin_count;
                     i < (n + 1) * (new_count / N) + left_marin_count; i++) {
                    if (igh[i] < 90 + delta && igh[i] > 90 - delta) {
                        count++;
                    }
                }
                if (count / (new_count / (double)N) > 1 - rho) {
                    igtest_pass_count_tmp++;
                }
                total_test_retio[axis][candidate] += count / (new_count / (double)N);
            }
            total_test_retio[axis][candidate] /= N;
            total_test_pass_count[axis][candidate] = igtest_pass_count_tmp;
            delete[] igh;
        }
        double max = 0;
        int max_pass_count = 0;
        int index = 0;

        for (int count = 0; count < NUM_CANDIDATE; count++) {
            // コマ内に侵入していなければカウントする
            if (invasion_test(axis, slc_final[axis][count].position, (int)slc_final[axis][count].pixels.size(), slc_final[axis][count].theta) == OK) {
                if (max_pass_count <= total_test_pass_count[axis][count]) {
                    max = total_test_retio[axis][count];
                    max_pass_count = total_test_pass_count[axis][count];
                    index = count;
                }
            }
        }
        axis_test_retio[axis] = max;
        axis_test_retio_index[axis] = index;
    }

    // 分割判定
    if (total_test_pass_count[0][axis_test_retio_index[0]] >= N - M && axis_test_retio[0] > axis_test_retio[1]) {
        slice_line = slc_final[0].at(axis_test_retio_index[0]);
#ifdef DEBUG
        cout << "position: " << slice_line.position << endl;
        cout << "slice tate" << endl;
#endif
    }
    else if (total_test_pass_count[1][axis_test_retio_index[1]] >= N - M && axis_test_retio[0] < axis_test_retio[1]) {
        slice_line = slc_final[1].at(axis_test_retio_index[1]);
#ifdef DEBUG
        cout << "position: " << slice_line.position << endl;
        cout << "slice yoko" << endl;
#endif
    }

    slc[0].clear();
    slc[1].clear();
}

void FrameSeparation::calculate_ig()
{
    // sobel filter
    IplImage* sobel_x = cvCloneImage(proc_img);
    IplImage* sobel_y = cvCloneImage(proc_img);
    cvSobel(src, sobel_x, 1, 0, 3);
    cvConvertScaleAbs(sobel_x, sobel_x, 1, 0);
    cvSobel(src, sobel_y, 0, 1, 3);
    cvConvertScaleAbs(sobel_y, sobel_y, 1, 0);
    cvAddWeighted(sobel_x, 0.5, sobel_y, 0.5, 0, proc_img);
    cvSmooth(proc_img, proc_img, CV_GAUSSIAN, 3);
    cvReleaseImage(&sobel_x);
    cvReleaseImage(&sobel_y);

#ifdef DEBUG
    cvNamedWindow("[ calculate_ig ] sobel image", cv::WINDOW_AUTOSIZE);
    cvShowImage("[ calculate_ig ] sobel image", proc_img);
    cv::waitKey(0);
#endif

    // width x height, 3chの行列
    this->ig_mat = Mat(CvSize(this->src->width, this->src->height), CV_MAKE_TYPE(CV_8U, 3));

    for (int y = (BLOCK_SIZE - 1) / 2; y < this->src->height - (BLOCK_SIZE - 1) / 2; y++) {
        for (int x = (BLOCK_SIZE - 1) / 2; x < this->src->width - (BLOCK_SIZE - 1) / 2; x++) {
            double fx = 0, fy = 0, magnitude = 0, direction = 0;
            fx = (int)(uchar) this->proc_img->imageData[this->proc_img->widthStep * y + (x + 1) * this->proc_img->nChannels] - (int)(uchar) this->proc_img->imageData[this->proc_img->widthStep * y + (x - 1) * this->proc_img->nChannels];
            fy = (int)(uchar) this->proc_img->imageData[this->proc_img->widthStep * (y + 1) + x * this->proc_img->nChannels] - (int)(uchar) this->proc_img->imageData[this->proc_img->widthStep * (y - 1) + x * this->proc_img->nChannels];
            magnitude = sqrt(fx * fx + fy * fy);
            direction = atan(fy / (fx + 0.01));
            direction = (direction < 0 ? direction + CV_PI : direction) * 180.0 / CV_PI;

            ig_mat.data[ig_mat.step * y + (x * ig_mat.channels())] = (uchar) this->src->imageData[this->src->widthStep * y + x * this->src->nChannels];
            ig_mat.at<Vec3b>(y, x)[1] = (uchar)direction;
        }
    }
#ifdef DEBUG
    vector<Mat> planes;
    split(ig_mat, planes);
    namedWindow("[ calculate_ig ] ig_mat[1]", cv::WINDOW_AUTOSIZE);
    imshow("[ calculate_ig ] ig_mat[1]", planes.at(1));
#endif
}

bool FrameSeparation::sl_exists()
{
    if (slice_line.position == 0) {
        return false;
    }
    if (slice_line.position <= (BLOCK_SIZE - 1) / 2) {
        return false;
    }
    int length = slice_line.is_horizontal ? src->height : src->width;
    int ep = slice_line.is_horizontal ? slice_line.pixels.at(slice_line.pixels.size() - 1).y : slice_line.pixels.at(slice_line.pixels.size() - 1).x;
    if (slice_line.theta != 90) {
        if (abs(length - slice_line.position) <= 20 || abs(length - ep) <= 20) {
            return false;
        }
    }
    return true;
}

int FrameSeparation::separation()
{
    rel_slice_point = rel_original_point;
    rel_remained_point = rel_original_point;

    vector<PixPoint>* pixels = &slice_line.pixels;
    // x軸方向に分割した場合
    if (slice_line.is_horizontal) {
        CvSize slice_size, remained_size;
        if (slice_line.theta == 90) {
            slice_size = CvSize(src->width, slice_line.position);
            remained_size = CvSize(src->width, src->height - slice_line.position);
        }
        else {
            if (pixels->at(pixels->size() - 1).y < slice_line.position) {
                slice_size = CvSize(pixels->at(pixels->size() - 1).x + 1, slice_line.position + 1);
            }
            else {
                slice_size = CvSize(pixels->at(pixels->size() - 1).x + 1,
                    pixels->at(pixels->size() - 1).y + 1);
            }
            if (pixels->at(pixels->size() - 1).y < slice_line.position) {
                remained_size = CvSize(
                    src->width + 1, src->height - pixels->at(pixels->size() - 1).y + 1);
            }
            else {
                remained_size = CvSize(src->width + 1, src->height - slice_line.position + 1);
            }
        }
        slice_src = cvCreateImage(slice_size, src->depth, src->nChannels);
        cvSet(slice_src, cvScalarAll(255), 0);
        remained_src = cvCreateImage(remained_size, src->depth, src->nChannels);
        cvSet(remained_src, cvScalarAll(255), 0);

        rel_remained_point.y += slice_line.position;

        int h, w, c;
        if (src->width == pixels->size() && slice_line.theta == 90) {
            for (w = 0; w < src->width; w++) {
                for (h = 0; h < src->height; h++) {
                    if (h < slice_line.position) {
                        for (c = 0; c < src->nChannels; c++) {
                            slice_src->imageData[slice_src->widthStep * h + w * slice_src->nChannels + c] = src->imageData[src->widthStep * h + w * src->nChannels + c];
                        }
                    }
                    else {
                        for (c = 0; c < src->nChannels; c++) {
                            remained_src->imageData[remained_src->widthStep * (h - slice_line.position) + w * remained_src->nChannels + c] = src->imageData[src->widthStep * h + w * src->nChannels + c];
                        }
                    }
                }
            }
        }
        else if (src->width == pixels->size() && slice_line.theta > 90) {
            for (w = 0; w < src->width; w++) {
                for (h = 0; h < src->height; h++) {
                    if (h < pixels->at(w).y) {
                        for (c = 0; c < src->nChannels; c++) {
                            slice_src->imageData[slice_src->widthStep * h + w * slice_src->nChannels + c] = src->imageData[src->widthStep * h + w * src->nChannels + c];
                        }
                    }
                    else {
                        for (c = 0; c < src->nChannels; c++) {
                            remained_src->imageData[remained_src->widthStep * (h - pixels->at(pixels->size() - 1).y) + w * slice_src->nChannels + c] = src->imageData[src->widthStep * h + w * src->nChannels + c];
                        }
                    }
                }
            }
        }
        else if (src->width == pixels->size() && slice_line.theta <= 90) {
            for (w = 0; w < src->width; w++) {
                for (h = 0; h < src->height; h++) {
                    if (h < pixels->at(w).y) {
                        for (c = 0; c < src->nChannels; c++) {
                            slice_src->imageData[slice_src->widthStep * h + w * slice_src->nChannels + c] = src->imageData[src->widthStep * h + w * src->nChannels + c];
                        }
                    }
                    else {
                        for (c = 0; c < src->nChannels; c++) {
                            remained_src->imageData[remained_src->widthStep * (h - slice_line.position) + w * slice_src->nChannels + c] = src->imageData[src->widthStep * h + w * src->nChannels + c];
                        }
                    }
                }
            }
        }
        else if (src->width > pixels->size() && slice_line.theta > 90) {
            for (w = 0; w < pixels->size(); w++) {
                for (h = 0; h < src->height; h++) {
                    if (h < pixels->at(w).y) {
                        for (c = 0; c < src->nChannels; c++) {
                            slice_src->imageData[slice_src->widthStep * h + w * slice_src->nChannels + c] = src->imageData[src->widthStep * h + w * src->nChannels + c];
                        }
                    }
                    else {
                        for (c = 0; c < src->nChannels; c++) {
                            remained_src->imageData[remained_src->widthStep * h + w * remained_src->nChannels + c] = src->imageData[src->widthStep * h + w * src->nChannels + c];
                        }
                    }
                }
            }
            for (w = (int)pixels->size(); w < src->width; w++) {
                for (h = 0; h < src->height; h++) {
                    for (c = 0; c < src->nChannels; c++) {
                        remained_src->imageData[remained_src->widthStep * h + w * remained_src->nChannels + c] = src->imageData[src->widthStep * h + w * src->nChannels + c];
                    }
                }
            }
        }
        else {
            for (w = 0; w < pixels->size(); w++) {
                for (h = 0; h < src->height; h++) {
                    if (h < pixels->at(w).y) {
                        for (c = 0; c < src->nChannels; c++) {
                            slice_src->imageData[slice_src->widthStep * h + w * slice_src->nChannels + c] = src->imageData[src->widthStep * h + w * src->nChannels + c];
                        }
                    }
                    else {
                        for (c = 0; c < src->nChannels; c++) {
                            remained_src->imageData[remained_src->widthStep * (h - slice_line.position) + w * remained_src->nChannels + c] = src->imageData[src->widthStep * h + w * src->nChannels + c];
                        }
                    }
                }
            }
            for (w = (int)pixels->size(); w < src->width; w++) {
                for (h = 0; h < src->height; h++) {
                    for (c = 0; c < src->nChannels; c++) {
                        slice_src->imageData[slice_src->widthStep * h + w * slice_src->nChannels + c] = src->imageData[src->widthStep * h + w * src->nChannels + c];
                    }
                }
            }
        }
    }
    // y軸方向に分割した場合
    else {
        CvSize slice_size, remained_size;
        if (slice_line.theta == 90) {
            slice_size = CvSize(slice_line.position, src->height);
            remained_size = CvSize(src->width - slice_line.position, src->height);
        }
        else {
            if (pixels->at(pixels->size() - 1).x < slice_line.position) {
                slice_size = CvSize(slice_line.position + 1, pixels->at(pixels->size() - 1).y + 1);
            }
            else {
                slice_size = CvSize(pixels->at(pixels->size() - 1).x + 1, src->height); // pixels->at(pixels->size() - 1).y+1);
            }
            if (pixels->at(pixels->size() - 1).x < slice_line.position) {
                remained_size = CvSize(src->width - pixels->at(pixels->size() - 1).x + 1, src->height);
            }
            else {
                remained_size = CvSize(src->width - slice_line.position, pixels->at(pixels->size() - 1).y + 1);
            }
        }
        slice_src = cvCreateImage(slice_size, src->depth, src->nChannels);
        cvSet(slice_src, cvScalarAll(255), 0);
        remained_src = cvCreateImage(remained_size, src->depth, src->nChannels);
        cvSet(remained_src, cvScalarAll(255), 0);

        rel_remained_point.x += slice_line.position;

        int h, w, c;
        if (src->height == pixels->size() && slice_line.theta == 90) {
            for (h = 0; h < src->height; h++) {
                for (w = 0; w < src->width; w++) {
                    if (w < slice_line.position) {
                        for (c = 0; c < src->nChannels; c++) {
                            slice_src->imageData[slice_src->widthStep * h + w * slice_src->nChannels + c] = src->imageData[src->widthStep * h + w * src->nChannels + c];
                        }
                    }
                    else {
                        for (c = 0; c < src->nChannels; c++) {
                            remained_src->imageData[remained_src->widthStep * h + (w - slice_line.position) * remained_src->nChannels + c] = src->imageData[src->widthStep * h + w * src->nChannels + c];
                        }
                    }
                }
            }
        }
        else if (src->height == pixels->size() && slice_line.theta > 90) {
            for (h = 0; h < src->height; h++) {
                for (w = 0; w < src->width; w++) {
                    if (w < pixels->at(h).x) {
                        for (c = 0; c < src->nChannels; c++) {
                            slice_src->imageData[slice_src->widthStep * h + w * slice_src->nChannels + c] = src->imageData[src->widthStep * h + w * src->nChannels + c];
                        }
                    }
                    else {
                        for (c = 0; c < src->nChannels; c++) {
                            remained_src->imageData[remained_src->widthStep * h + (w - pixels->at(pixels->size() - 1).x) * slice_src->nChannels + c] = src->imageData[src->widthStep * h + w * src->nChannels + c];
                        }
                    }
                }
            }
        }
        else if (src->height == pixels->size() && slice_line.theta < 90) {
            for (h = 0; h < src->height; h++) {
                for (w = 0; w < src->width; w++) {
                    if (w < pixels->at(h).x) {
                        for (c = 0; c < src->nChannels; c++) {
                            slice_src->imageData[slice_src->widthStep * h + w * slice_src->nChannels + c] = src->imageData[src->widthStep * h + w * src->nChannels + c];
                        }
                    }
                    else {
                        for (c = 0; c < src->nChannels; c++) {
                            remained_src->imageData[remained_src->widthStep * h + (w - slice_line.position) * slice_src->nChannels + c] = src->imageData[src->widthStep * h + w * src->nChannels + c];
                        }
                    }
                }
            }
        }
        else if (src->height > pixels->size() && slice_line.theta > 90) {
            for (h = 0; h < pixels->size(); h++) {
                for (w = 0; w < src->width; w++) {
                    if (w < pixels->at(h).x) {
                        for (c = 0; c < src->nChannels; c++) {
                            slice_src->imageData[slice_src->widthStep * h + w * slice_src->nChannels + c] = src->imageData[src->widthStep * h + w * src->nChannels + c];
                        }
                    }
                    else {
                        for (c = 0; c < src->nChannels; c++) {
                            remained_src->imageData[remained_src->widthStep * h + w * remained_src->nChannels + c] = src->imageData[src->widthStep * h + w * src->nChannels + c];
                        }
                    }
                }
            }
            for (h = (int)pixels->size(); h < src->height; h++) {
                for (w = 0; w < src->width; w++) {
                    for (c = 0; c < src->nChannels; c++) {
                        remained_src->imageData[remained_src->widthStep * h + w * remained_src->nChannels + c] = src->imageData[src->widthStep * h + w * src->nChannels + c];
                    }
                }
            }
        }
        else {
            for (h = 0; h < pixels->size(); h++) {
                for (w = 0; w < src->width; w++) {
                    if (w < pixels->at(h).x) {
                        for (c = 0; c < src->nChannels; c++) {
                            slice_src->imageData[slice_src->widthStep * h + w * slice_src->nChannels + c] = src->imageData[src->widthStep * h + w * src->nChannels + c];
                        }
                    }
                    else {
                        for (c = 0; c < src->nChannels; c++) {
                            remained_src->imageData[remained_src->widthStep * h + (w - slice_line.position) * remained_src->nChannels + c] = src->imageData[src->widthStep * h + w * src->nChannels + c];
                        }
                    }
                }
            }
            for (h = (int)pixels->size(); h < src->height; h++) {
                for (w = 0; w < src->width; w++) {
                    for (c = 0; c < src->nChannels; c++) {
                        remained_src->imageData[remained_src->widthStep * h + w * remained_src->nChannels + c] = src->imageData[src->widthStep * h + w * src->nChannels + c];
                    }
                }
            }
        }
    }
#ifdef DEBUG
    cvDestroyWindow("[ separation ] slice_src");
    cvShowImage("[ separation ] slice_src", slice_src);
    cv::waitKey(0);
    cvDestroyWindow("[ separation ] remained_src");
    cvShowImage("[ separation ] remained_src", remained_src);
    cv::waitKey(0);
#endif

    // 保存しないで良い余白等はfalseを返す
    if (is_blank(slice_src)) {
        return DROP_SLICE_SRC;
    }
    if (is_blank(remained_src)) {
        return DROP_REMAINED_SRC;
    }

    double threshold = 0.02;
    if (slice_src->width * slice_src->height < original_size * threshold) {
        return DROP_SLICE_SRC;
    }
    if (remained_src->width * remained_src->height < original_size * threshold) {
        return DROP_REMAINED_SRC;
    }

    return OK;
}

bool FrameSeparation::is_blank(IplImage* image)
{
    int count = 0;
    int area = image->width * image->height;
    for (int h = 0; h < image->height; h++) {
        for (int w = 0; w < image->width; w++) {
            if ((int)(unsigned char)image->imageData[image->widthStep * h + w * image->nChannels] < 200) {
                count++;
            }
        }
    }
// 黒の割合
#ifdef DEBUG
    cout << "black area rate: " << count / (double)area << endl;
#endif
    if (count / (double)area < 0.05) {
        return true;
    }
    if (src->width * src->height < original_size * 0.05) {
        return true;
    }
    return false;
}

// is_horizontal : 分割の切り口が水平かどうか
void FrameSeparation::calculate_slc(bool is_horizontal)
{
    int position = is_horizontal ? proc_img->height : proc_img->width;
    int direction = is_horizontal ? proc_img->width : proc_img->height;

#ifdef DEBUG
    double* ig_line = new double[position];
#endif
    vector<PixPoint>* line;
    slc[is_horizontal].resize(position);

    for (int p = (BLOCK_SIZE - 1) / 2; p < position - (BLOCK_SIZE - 1) / 2; p++) {
        double ig = 0;
        int num_zero = 0;
        line = &slc[is_horizontal].at(p).pixels;
        SL* sl = &slc[is_horizontal].at(p);

        detect_pixels(is_horizontal, p, direction, 90, line);
        for (int d = (BLOCK_SIZE - 1) / 2; d < direction - (BLOCK_SIZE - 1) / 2;
             d++) {
            if (is_horizontal) {
                // 90度が一番高くなるように正規化する
                if (ig_mat.at<Vec3b>(line->at(d).y, line->at(d).x)[1] >= 90) {
                    ig += 180 - ig_mat.at<Vec3b>(line->at(d).y, line->at(d).x)[1];
                }
                else {
                    ig += ig_mat.at<Vec3b>(line->at(d).y, line->at(d).x)[1];
                }
            }
            else {

                // 0, 180度が一番高くなるように正規化する
                if (ig_mat.at<Vec3b>(line->at(d).y, line->at(d).x)[1] >= 90) {
                    ig += -90 + ig_mat.at<Vec3b>(line->at(d).y, line->at(d).x)[1];
                }
                else {
                    ig += 90 - ig_mat.at<Vec3b>(line->at(d).y, line->at(d).x)[1];
                }
            }
        }
#ifdef DEBUG
        ig_line[p] = ig / line->size();
#endif
        sl->is_horizontal = is_horizontal;
        sl->position = p;
        sl->theta = 90;
        sl->ig = num_zero > line->size() * 0.9 ? 0 : ig / line->size();
    }
}

void FrameSeparation::calculate_wpr(bool is_horizontal)
{
    int position = is_horizontal ? proc_img->height : proc_img->width;
    int direction = is_horizontal ? proc_img->width : proc_img->height;

#ifdef DEBUG
    double* wpr = new double[position];
#endif
    Mat wpr_line = Mat::zeros(direction, 1, CV_8U);

    for (int p = (BLOCK_SIZE - 1) / 2; p < position - (BLOCK_SIZE - 1) / 2; p++) {
        // wprを求める
        int max = 0, tmp = 0;
        for (int d = 0; d < direction; d++) {
            if (is_horizontal) {
                wpr_line.at<uchar>(0, d) = (uchar)proc_img->imageData[proc_img->widthStep * p + d * proc_img->nChannels];
            }
            else {
                wpr_line.at<uchar>(0, d) = (uchar)proc_img->imageData[proc_img->widthStep * d + p * proc_img->nChannels];
            }
        }
        for (int i = 0; i < direction; i++) {
            if (wpr_line.at<uchar>(0, i) == 255) {
                tmp++;
            }
            else {
                max = (max > tmp ? max : tmp);
                tmp = 0;
            }
        }
        max = (max > tmp ? max : tmp);
#ifdef DEBUG
        wpr[p] = max / (double)direction;
#endif
        slc[is_horizontal].at(p).wpr = max / (double)direction;
    }
}

// is_horizontal : 分割の切り口が水平かどうか
void FrameSeparation::calculate_oblique_slc(bool is_horizontal)
{
    // 走査する対象の線の画素情報を格納する
    vector<PixPoint>* line;

    int position = is_horizontal ? proc_img->height : proc_img->width;
    int direction = is_horizontal ? proc_img->width : proc_img->height;

    // 水平・垂直方向のデータが残っているので初期化する
    slc[is_horizontal].clear();

    // 45度〜135度までの直線を走査する
    slc[is_horizontal].resize(position * (135 - 45));
#ifdef DEBUG
    double* ig_line = new double[position * (135 - 45)];
#endif
    SL* sl;
    int num_zero, theta, d;
    double ig;
    for (int p = (BLOCK_SIZE - 1) / 2; p < position - (BLOCK_SIZE - 1) / 2; p++) {
        int count = 0;
        for (theta = 45; theta < 135; theta++) {
            num_zero = 0;
            sl = &slc[is_horizontal].at(p * (135 - 45) + count);
            line = &sl->pixels;
            try {
                detect_pixels(is_horizontal, p, direction, theta, line);
            }
            catch (Response err) {
                if (err == INVALID_DIGREES) {
#ifdef DEBUG
                    cout << "invalid digree" << endl;
#endif
                }
                count++;
                continue;
            }

            if (line->size() < direction * 0.8) {
                count++;
                continue;
            }
            ig = 0.0;
            for (d = 0; d < line->size(); d++) {
                if (is_horizontal) {
                    if (ig_mat.data[line->at(d).y * ig_mat.step + line->at(d).x * ig_mat.elemSize() + 1] >= 90) {
                        ig += 180 - ig_mat.data[line->at(d).y * ig_mat.step + line->at(d).x * ig_mat.elemSize() + 1];
                    }
                    else {
                        ig += ig_mat.data[line->at(d).y * ig_mat.step + line->at(d).x * ig_mat.elemSize() + 1];
                    }
                }
                else {
                    if (ig_mat.data[line->at(d).y * ig_mat.step + line->at(d).x * ig_mat.elemSize() + 1] == 0) {
                        num_zero++;
                    }
                    // 0, 180度が一番高くなるように正規化する
                    if (ig_mat.data[line->at(d).y * ig_mat.step + line->at(d).x * ig_mat.elemSize() + 1] >= 90) {
                        ig += -90 + ig_mat.data[line->at(d).y * ig_mat.step + line->at(d).x * ig_mat.elemSize() + 1];
                    }
                    else {
                        ig += 90 - ig_mat.data[line->at(d).y * ig_mat.step + line->at(d).x * ig_mat.elemSize() + 1];
                    }
                }
            }
#ifdef DEBUG
            ig_line[p * (135 - 45) + count] = ig / line->size();
#endif
            sl->is_horizontal = is_horizontal;
            sl->position = p;
            sl->theta = theta;
            sl->ig = num_zero > line->size() * 0.9 ? 0 : ig / line->size();
            count++;
        }
    }
}

void FrameSeparation::calculate_oblique_wpr(bool is_horizontal)
{
    // 走査する対象の線の画素情報を格納する
    vector<PixPoint>* line;

    int position = is_horizontal ? proc_img->height : proc_img->width;
    int direction = is_horizontal ? proc_img->width : proc_img->height;

#ifdef DEBUG
    double* wpr = new double[position * (135 - 45)];
#endif
    int max = 0, tmp = 0, count = 0;
    for (int p = (BLOCK_SIZE - 1) / 2; p < position - (BLOCK_SIZE - 1) / 2; p++) {
        count = 0;
        for (int theta = 45; theta < 135; theta++) {
            // wprを求める
            max = 0, tmp = 0;
            line = &slc[is_horizontal].at((p) * (135 - 45) + count).pixels;
            if (line->size() < direction * 0.5) {
                count++;
                continue;
            }
            for (int i = 0; i < line->size(); i++) {
                if ((uchar)bin_img->imageData[bin_img->widthStep * line->at(i).y + line->at(i).x * bin_img->nChannels] == 255) {
                    tmp++;
                }
                else {
                    max = max > tmp ? max : tmp;
                    tmp = 0;
                }
            }
            max = max > tmp ? max : tmp;
#ifdef DEBUG
            wpr[(p) * (135 - 45) + count] = max / (double)line->size();
#endif
            slc[is_horizontal].at((p) * (135 - 45) + count).wpr = max / (double)line->size();

            count++;
        }
    }
}
