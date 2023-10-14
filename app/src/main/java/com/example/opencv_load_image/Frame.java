package com.example.opencv_load_image;

import org.opencv.core.Mat;
import org.opencv.core.MatOfPoint;
import org.opencv.core.MatOfPoint2f;
import org.opencv.core.Point;
import org.opencv.core.Size;
import org.opencv.imgproc.Imgproc;

import java.util.ArrayList;
import java.util.List;

public class Frame {
//                v_idx = -1
//                tot = 0
//                div = 0
//                impulses = []
//                while v_idx < len(vals)-1:
//                v_idx += 1
//                if vals[v_idx] < (.3 * max(vals)):
//                continue
//                while vals[v_idx] > (.3 * max(vals)):
//                tot += v_idx * vals[v_idx]
//                div += vals[v_idx]
//                v_idx += 1
//                try:
//                impulses.append(float(tot) / div)
//                except ZeroDivisionError:
//                pass
//                        tot = 0
//                div = 0


//                Mat drawing = Mat.zeros(mat.size(), CvType.CV_8UC3);
//                List<MatOfPoint> contoursPolyList = new ArrayList<>(contoursPoly.length);
//
//                for (MatOfPoint2f poly : contoursPoly) {
//                    contoursPolyList.add(new MatOfPoint(poly.toArray()));
//                }
//
//                Scalar color = new Scalar(rng.nextInt(256), rng.nextInt(256), rng.nextInt(256));
//                for (int i = 0; i < contours.size(); i++) {
//                    Imgproc.drawContours(drawing, contoursPolyList, i, color);
//                }
//

//                Imgproc.circle(drawing, centers[index], (int) radius[index][0], color, 2);
//
//                Utils.matToBitmap(drawing, bitmap);
    public Mat mat;
    public double threshold;

    public Frame(Mat mat, double threshold){
        this.mat = mat;
        this.threshold = threshold;
    }

    public static Mat find_contour(List<MatOfPoint> contours, Mat mat){
        Mat hierarchy = new Mat();
        Mat mat_copy;
        Imgproc.cvtColor(mat, mat, Imgproc.COLOR_RGB2GRAY);
        mat_copy = mat.clone();
        Imgproc.blur(mat, mat, new Size(5, 5));
        Imgproc.threshold(mat, mat ,0, 255, Imgproc.THRESH_BINARY + Imgproc.THRESH_OTSU);
        Imgproc.findContours(mat, contours, hierarchy, Imgproc.RETR_LIST, Imgproc.CHAIN_APPROX_SIMPLE);
        hierarchy.release();
        return mat_copy;
    }

    public static double[] radius_and_center(List<MatOfPoint> contours){
        MatOfPoint2f[] contoursPoly = new MatOfPoint2f[contours.size()];
        Point[] centers = new Point[contours.size()];
        float[][] radius = new float[contours.size()][1];
        for (int i = 0; i < contours.size(); i++) {
            contoursPoly[i] = new MatOfPoint2f();
            Imgproc.approxPolyDP(new MatOfPoint2f(contours.get(i).toArray()), contoursPoly[i], 3, true);
            centers[i] = new Point();
            Imgproc.minEnclosingCircle(contoursPoly[i], centers[i], radius[i]);
        }

        float max_r = 0;
        int index = 0;
        for (int i = 0; i < radius.length; i++) {
            if (i == 0)
                max_r = radius[i][0];
            else {
                if (radius[i][0] > max_r) {
                    max_r = radius[i][0];
                    index = i;
                }
            }
        }

        double[] r_n_c = {(double) max_r, centers[index].x, centers[index].y};
        return r_n_c;
    }

    public static int[] frame_quantify(Mat edge, int[] range, float radius){

        int[] columnSums = new int[2 * (int) radius + 40];
        for (int col = range[0]; col < range[1]; col++) {
            int sum = 0;
            for (int row = range[2]; row < range[3]; row++) {
                sum += (int) edge.get(row, col)[0]; // Access the pixel value at (row, col)
            }
            columnSums[col- range[0]] = sum;
        }

        return columnSums;
    }

    public static List<Float> getPulses(int[] columnSums, int start_x) {
        int column_index = -1;
        int weighted_sum = 0;
        int sum_weight = 0;
        List<Float> pulse = new ArrayList<>();
        int max = 0;
        for (int value : columnSums) {
            if (value > max)
                max = value;
        }
        while (column_index < columnSums.length - 1) {
            column_index++;
            if ((float) columnSums[column_index] > 0.2 * max) {
                while ((float) columnSums[column_index] > 0.2 * max) {
                    weighted_sum += (column_index + start_x) * columnSums[column_index];
                    sum_weight += columnSums[column_index];
                    column_index++;
                }
                pulse.add((float) weighted_sum / (float) sum_weight);
                weighted_sum = 0;
                sum_weight = 0;
            }
        }
        return pulse;
    }

    public static void filterPules(List<Float> pulse, Mat mat, double threshold, int center_y, float radius){
        boolean last = false;
        boolean cur = false;
        List<Integer> index = new ArrayList<>();
        int count = 0;
        for(int i = 1; i < pulse.size(); i++){
            float check = pulse.get(i - 1) + (pulse.get(i) - pulse.get(i - 1)) / 2;

            if ((double) mat.get(center_y, (int) check)[0] > threshold)
                count++;
            if ((double) mat.get(center_y - (int) (radius / 2), (int) check)[0] > threshold)
                count++;
            if ((double) mat.get(center_y + (int) (radius / 2), (int) check)[0] > threshold)
                count++;

            if (count >= 2)
                cur = true;
            if (last && cur)
                index.add(i);
            last = cur;
            cur = false;
            count = 0;
        }
        for (int i = 0; i < index.size(); i++){
            pulse.remove(index.get(i) - i);
        }
    }
    public static Frame frame_process(Mat mat_copy, float radius){
        float factor = 1 / radius;
        int blur_width = (int) (factor * radius + 1);
        Mat inter = new Mat();
        Imgproc.blur(mat_copy, inter, new Size(blur_width, radius));
        double threshold = Imgproc.threshold(inter, inter ,0, 255,
                Imgproc.THRESH_BINARY + Imgproc.THRESH_OTSU);

        return new Frame(inter, threshold);
    }

    public static Mat draw(Mat mat, List<Float> pulse, int center_y){
        int numRows = mat.rows();
        double[] a = {200};
        for (int i = 0; i < pulse.size(); i++){
            for (int y = 0; y < numRows; y++){
                int x = 0;
                if(pulse.get(i) % ((float) ((int) pulse.get(1).floatValue())) >= 0.5)
                    x = (int) pulse.get(i).floatValue() + 1;
                else
                    x = (int) pulse.get(i).floatValue();

                mat.put(y, x, a);
            }
        }
        for (int i = 0; i < mat.cols(); i++){
            mat.put(center_y, i, a);
        }
        return mat;
    }

    public static float calc_freq(List<Float> pulse){
        List<Float> pulse_copy_1  = new ArrayList<>();
        List<Float> pulse_copy_2  = new ArrayList<>();
        for (float nums: pulse){
            pulse_copy_1.add(nums);
            pulse_copy_2.add(nums);
        }
        pulse_copy_1.remove(0);
        pulse_copy_2.remove(pulse_copy_2.size() - 1);
        List<Float> diff = new ArrayList<>();
        float sum = 0;
        for(int i = 0; i < pulse_copy_2.size(); i++){
            sum += pulse_copy_1.get(i) - pulse_copy_2.get(i);
        }
        float interval = sum / pulse_copy_2.size();
//        float period = interval * 0.0000319f;
        float period = interval * 0.00001858f;
        float freq = 1 / period;
        return freq;
    }

    public static int[] light_range(int center_x, int center_y, float radius){
        int start_x = center_x - (int) radius - 20;
        int end_x = center_x + (int) radius + 20;
        int start_y = center_y - (int) radius - 20;
        int end_y = center_y + (int) radius + 20;
        return new int[]{start_x, end_x, start_y, end_y};
    }
}
