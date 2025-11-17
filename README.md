# Mô phỏng so sánh Bộ lọc Wiener và Bộ lọc LMS

Đây là dự án MATLAB mô phỏng và so sánh hiệu quả của bộ lọc Wiener và bộ lọc LMS trong việc xử lý tín hiệu.

##  Yêu cầu Cài đặt

* MATLAB 
##  Hướng dẫn Chạy
Bước 1 : Tải folder data, cho folder data vào cùng workplace với file mix_data.m , sau khi chạy file mix_data.m sẽ cho ra folder mixed chứa âm thanh trộn nhiễu ở folder clean và folder noise.

Bước 2 : Cho folder âm thanh mixed, clean và noise cùng trong workplace mới, cho chạy file lms_filter và weiner_filter. Khi chạy file lms_filter.m sẽ cho folder âm thanh đã qua bộ lọc lms với tên folder lms_filtered_fix.
Tương tự với file weiner.m sẽ cho folder âm thanh đã qua bộ lọc Weiner với tên folder weiner_filtered_fix.

Bước 3 : Cho folder lms_filtered_fix, folder weiner_filtered_fix vào cùng workplace với 2 file weiner_sosanh.m và lms_sosanh.m, khi chạy file weiner_sosanh.m sẽ cho ra biểu đồ dạng sóng - phổ với 3 trường hợp NOISY - FILTERED -
CLEAN so sánh trước khi khử nhiễu, sau khi khử nhiễu , dữ liệu gốc cùng với file csv weiner_compare_metrics so sánh delta SNR và tỷ lệ hiệu quả khử nhiễu ở các trường hợp SNR đầu vào khác nhau.

Bước 4 : Chạy file sosanh2pp.m , cho ra kết quả sánh SNR và tỷ lệ hiệu quả khử nhiễu giữa bộ lọc Wiener và LMS, đồng thời xuất ra biểu đồ so sánh 2 bộ lọc.

Bước 5 : Chạy file kiemtraμ.m so sánh các μ khác nhau trong bộ lọc LMS, xuất ra được biểu đồ LMS ở các μ khác nhau. 
