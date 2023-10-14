package com.example.opencv_load_image;

import android.content.Intent;
import android.graphics.Bitmap;
import android.graphics.PointF;
import android.media.MediaMetadataRetriever;
import android.net.Uri;
import android.os.Build;
import android.os.Bundle;
import android.os.Handler;
import android.provider.MediaStore;
import android.util.Log;
import android.view.MotionEvent;
import android.view.View;
import android.view.Window;
import android.view.WindowManager;
import android.widget.Button;
import android.widget.ImageView;
import android.widget.TextView;

import androidx.annotation.Nullable;
import androidx.appcompat.app.AppCompatActivity;

import org.opencv.android.BaseLoaderCallback;
import org.opencv.android.LoaderCallbackInterface;
import org.opencv.android.OpenCVLoader;
import org.opencv.android.Utils;
import org.opencv.core.Core;
import org.opencv.core.CvType;
import org.opencv.core.Mat;
import org.opencv.core.MatOfPoint;
import org.opencv.core.MatOfPoint2f;
import org.opencv.core.Point;
import org.opencv.core.Rect;
import org.opencv.core.Scalar;
import org.opencv.core.Size;
import org.opencv.imgproc.Imgproc;

import java.io.IOException;
import java.net.URI;
import java.util.ArrayList;
import java.util.List;
import java.util.Random;

public class MainActivity extends AppCompatActivity {
    Button select;
    Button video;
    ImageView imageView;
    TextView textView;
    Bitmap bitmap;

    MediaMetadataRetriever mediaMetadataRetriever;
    final Handler handler = new Handler();
    float[] lastEvent = null;
    float d = 0f;
    float newRot = 0f;
    private boolean isZoomAndRotate;
    private boolean isOutSide;
    private static final int NONE = 0;
    private static final int DRAG = 1;
    private static final int ZOOM = 2;
    private int mode = NONE;
    private PointF start = new PointF();
    private PointF mid = new PointF();
    float oldDist = 1f;
    private float xCoOrdinate, yCoOrdinate;




    private BaseLoaderCallback mLoaderCallback = new BaseLoaderCallback(this) {
        @Override
        public void onManagerConnected(int status) {
            switch (status) {
                case LoaderCallbackInterface.SUCCESS:
                {
                    Log.i("OpenCV", "OpenCV loaded successfully");

                } break;
                default:
                {
                    super.onManagerConnected(status);
                } break;
            }
        }
    };

    @Override
    protected void onCreate(Bundle savedInstanceState){
        super.onCreate(savedInstanceState);
        setContentView(R.layout.layout);

        mediaMetadataRetriever = new MediaMetadataRetriever();

        select = findViewById(R.id.select);
        video = findViewById(R.id.video);
        imageView = findViewById(R.id.image);
        textView = findViewById(R.id.freq);

        imageView.setOnTouchListener(new View.OnTouchListener() {
            @Override
            public boolean onTouch(View v, MotionEvent event) {
                ImageView view = (ImageView) v;
                view.bringToFront();
                viewTransformation(view, event);
                return true;
            }
        });

        select.setOnClickListener(new View.OnClickListener() {
            @Override
            public void onClick(View v) {
                Intent intent = new Intent(Intent.ACTION_GET_CONTENT);
                intent.setType("image/*");
                startActivityForResult(intent, 100);
            }
        });

        video.setOnClickListener(new View.OnClickListener() {
            @Override
            public void onClick(View v) {
                Intent intent = new Intent(Intent.ACTION_GET_CONTENT, MediaStore.Images.Media.EXTERNAL_CONTENT_URI);
                intent.setType("video/*");
                startActivityForResult(intent, 10);
            }
        });
    }

    @Override
    protected void onActivityResult(int requestCode, int resultCode, @Nullable Intent data) {
        super.onActivityResult(requestCode, resultCode, data);

        if(requestCode == 10 && data != null){
            new Thread(new Runnable() {
                public void run() {
                    // Perform image processing here
                    Uri vidUri = data.getData();
                    MediaMetadataRetriever tRetriever = new MediaMetadataRetriever();
                    tRetriever.setDataSource(getBaseContext(), vidUri);
                    mediaMetadataRetriever = tRetriever;
                    String dur = mediaMetadataRetriever.extractMetadata(
                            MediaMetadataRetriever.METADATA_KEY_DURATION);
                    float duration = Float.parseFloat(dur);

                    float freq_sum = 0;
                    int step = 5;
                    int frame_num = (int) duration / 33 - 1;
                    int frame_index = 0;
                    int offset = 0;
                    String freqs = "";

                    Mat mat_copy;
                    Mat mat = new Mat();
                    Mat edge = new Mat();
                    Mat result;

                    while(frame_index < frame_num){
                        bitmap = mediaMetadataRetriever.getFrameAtIndex(frame_index);

                        Utils.bitmapToMat(bitmap, mat);

                        List<MatOfPoint> contours = new ArrayList<>();
                        // Find contour of the led
                        mat_copy = Frame.find_contour(contours, mat);
                        Utils.matToBitmap(mat_copy, bitmap);

                        // Calculate radius and center of the led
                        double[] radius_n_center = Frame.radius_and_center(contours);
                        int center_y = (int) radius_n_center[2];
                        int start = (int) (radius_n_center[1]);
                        int end = mat_copy.cols();
                        boolean black = true;

                        for (int i = start; i < end; i++){
                            double test = mat_copy.get(center_y, i)[0];
                            if (test > 50) {
                                black = false;
                                break;
                            }
                        }

                        if(!black) {
                            int[] range = Frame.light_range((int) radius_n_center[1], (int) radius_n_center[2], (float) radius_n_center[0]);
                            // Process image to get the edge of pulses
                            Frame frame = Frame.frame_process(mat_copy, (float) radius_n_center[0]);

                            Imgproc.Canny(frame.mat, edge, frame.threshold / 5, frame.threshold, 7);

                            // Convert image to quantifiable data
                            List<Float> pulse = Frame.getPulses(Frame.frame_quantify(edge, range, (float) radius_n_center[0]), range[0]);

                            Frame.filterPules(pulse, frame.mat, frame.threshold, (int) radius_n_center[2], (float) radius_n_center[0]);

                            // Calculate corresponding led frequency
                            float freq = Frame.calc_freq(pulse);
                            freq_sum += freq;
                            frame_index += step;
                            if(frame_index - offset == 30 || frame_index - offset == 60) {
                                freqs += Float.toString(freq_sum / 6) + "\n";
                                freq_sum = 0;
                            }

                            if (frame_index - offset >= 90) {
                                // Draw the edge of pulse
                                freqs += Float.toString(freq_sum / 6) + "\n";
                                result = Frame.draw(frame.mat, pulse, (int) radius_n_center[2]);
                                Utils.matToBitmap(result, bitmap);

                                break;
                            }
                        }else {
                            frame_index++;
                            offset++;
                        }
                    }

                    String finalFreqs = freqs;
                    runOnUiThread(new Runnable() {
                        public void run() {
                            // Update the UI with the processed image or results
                            imageView.setImageBitmap(bitmap);
                            textView.setText(finalFreqs);
                        }
                    });
                }
            }).start();

        }

        if(requestCode == 100 && data != null){
            Mat mat_copy;
            Mat result;
            Mat mat = new Mat();
            try {
                bitmap = MediaStore.Images.Media.getBitmap(this.getContentResolver(), data.getData());

                List<MatOfPoint> contours = new ArrayList<>();

                Utils.bitmapToMat(bitmap, mat);
                // Find contour of the led
                mat_copy = Frame.find_contour(contours, mat);

                // Calculate radius and center of the led
                double[] radius_n_center = Frame.radius_and_center(contours);
                int y = (int) radius_n_center[2];
                boolean black = true;
                for (int i = 0; i < mat_copy.cols(); i++){
                    double test = mat_copy.get(y, i)[0];
                    if (test != 0){
                        black = false;
                        break;
                    }
                }
                if (!black){
                    int[] range = Frame.light_range((int) radius_n_center[1], (int) radius_n_center[2], (float) radius_n_center[0]);

                    // Process image to get the edge of pulses
                    Frame struct =Frame.frame_process(mat_copy, (float) radius_n_center[0]);
                    Mat edge = new Mat();
                    Imgproc.Canny(struct.mat, edge, struct.threshold / 5, struct.threshold, 7);

                    // Convert image to quantifiable data
                    int[] columnSums = Frame.frame_quantify(edge, range, (float) radius_n_center[0]);
                    List<Float> pulse = Frame.getPulses(columnSums, range[0]);

                    Frame.filterPules(pulse, struct.mat, struct.threshold, (int) radius_n_center[2], (float) radius_n_center[0]);
                    // Calculate corresponding led frequency
                    float freq = Frame.calc_freq(pulse);

                    // Draw the edge of pulse
                    result = Frame.draw(struct.mat, pulse, (int) radius_n_center[2]);

                    Utils.matToBitmap(result, bitmap);
                    imageView.setImageBitmap(bitmap);
                    textView.setText(Float.toString(freq));
                }
            } catch (IOException e) {
                throw new RuntimeException(e);
            }

        }
    }
    public void onResume()
    {
        super.onResume();
        if (!OpenCVLoader.initDebug()) {
            Log.d("OpenCV", "Internal OpenCV library not found. Using OpenCV Manager for initialization");
            OpenCVLoader.initAsync(OpenCVLoader.OPENCV_VERSION_3_0_0, this, mLoaderCallback);
        } else {
            Log.d("OpenCV", "OpenCV library found inside package. Using it!");
            mLoaderCallback.onManagerConnected(LoaderCallbackInterface.SUCCESS);
        }
    }

    private void frame_process(){

    }
    private void viewTransformation(View view, MotionEvent event) {
        switch (event.getAction() & MotionEvent.ACTION_MASK) {
            case MotionEvent.ACTION_DOWN:
                xCoOrdinate = view.getX() - event.getRawX();
                yCoOrdinate = view.getY() - event.getRawY();

                start.set(event.getX(), event.getY());
                isOutSide = false;
                mode = DRAG;
                lastEvent = null;
                break;
            case MotionEvent.ACTION_POINTER_DOWN:
                oldDist = spacing(event);
                if (oldDist > 10f) {
                    midPoint(mid, event);
                    mode = ZOOM;
                }

                lastEvent = new float[4];
                lastEvent[0] = event.getX(0);
                lastEvent[1] = event.getX(1);
                lastEvent[2] = event.getY(0);
                lastEvent[3] = event.getY(1);
                d = rotation(event);
                break;
            case MotionEvent.ACTION_UP:
                isZoomAndRotate = false;
                if (mode == DRAG) {
                    float x = event.getX();
                    float y = event.getY();
                }
            case MotionEvent.ACTION_OUTSIDE:
                isOutSide = true;
                mode = NONE;
                lastEvent = null;
            case MotionEvent.ACTION_POINTER_UP:
                mode = NONE;
                lastEvent = null;
                break;
            case MotionEvent.ACTION_MOVE:
                if (!isOutSide) {
                    if (mode == DRAG) {
                        isZoomAndRotate = false;
                        view.animate().x(event.getRawX() + xCoOrdinate).y(event.getRawY() + yCoOrdinate).setDuration(0).start();
                    }
                    if (mode == ZOOM && event.getPointerCount() == 2) {
                        float newDist1 = spacing(event);
                        if (newDist1 > 10f) {
                            float scale = newDist1 / oldDist * view.getScaleX();
                            view.setScaleX(scale);
                            view.setScaleY(scale);
                        }
                        if (lastEvent != null) {
                            newRot = rotation(event);
                            view.setRotation((float) (view.getRotation() + (newRot - d)));
                        }
                    }
                }
                break;
        }
    }

    private float rotation(MotionEvent event) {
        double delta_x = (event.getX(0) - event.getX(1));
        double delta_y = (event.getY(0) - event.getY(1));
        double radians = Math.atan2(delta_y, delta_x);
        return (float) Math.toDegrees(radians);
    }

    private float spacing(MotionEvent event) {
        float x = event.getX(0) - event.getX(1);
        float y = event.getY(0) - event.getY(1);
        return (int) Math.sqrt(x * x + y * y);
    }

    private void midPoint(PointF point, MotionEvent event) {
        float x = event.getX(0) + event.getX(1);
        float y = event.getY(0) + event.getY(1);
        point.set(x / 2, y / 2);
    }

}
