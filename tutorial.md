# Tải mã nguồn:
```bash
git clone https://github.com/ailab-bio/GraphINVENT2.git
cd GraphINVENT2
```

# Install dependencies:
Lưu ý: Môi trường Conda sẽ có tên là graphinvent2
```bash
conda create -n graphinvent2 python=3.9 pytorch::pytorch torchvision torchaudio -c pytorch -y
conda activate graphinvent2
conda install -n graphinvent2 conda-forge::rdkit conda-forge::tqdm anaconda::h5py anaconda::scikit-learn matplotlib tensorboard -y
```

Error: Nếu gặp lỗi bên dưới khi chạy submit.py
```
ImportError: numpy.core.multiarray failed to import
AttributeError: _ARRAY_API not found
```

Bạn có thể thử cài đặt lại numpy với phiên bản thấp hơn 2.0

```bash
pip install "numpy<2.0"
```

# Tải dữ liệu:

Tải và giải nén dữ liệu từ:
```bash
cd /root/projects
wget https://zenodo.org/records/13897048/files/coconut_complete-10-2024.csv.zip?download=1
unzip <Tên tệp zip đã tải về>
```

Chuyển dữ liệu vào thư mục data:
```bash
mv <Tệp .csv đã giải nén> /root/projects/GraphINVENT2/data/
```

# Chạy mã nguồn:
## GDB13
Test thử với gdb13-debug với các thông số trong file submut_gdb13_debug.py. Chạy test thử với gdb13-debug:
```bash
python submit_gdb13_debug.py
``` 

## Coconut
Bước 1: Lọc dữ liệu từ file CSV, mục tiêu tạo đầu ra là tệp .smi mới với chỉ 2 cột: SMILES và Name. 
Thời gian chạy qui trình này phụ thuộc vào kích thước tệp đầu vào.
Đầu vào mặc định là tệp: data/pre-training/coconut/raw/coconut_complete-10-2024.csv
Tệp đầu ra sẽ được lưu trong thư mục data/pre-training/coconut/filter/coconut-complete-10-2024.smi
Nếu muốn đổi thành tệp khác, bạn vào mã nguồn và sửa lại tên tệp của thông số `input_file` và `output_file` trong file filter/filter_coconut_csv.py

```bash
python filter/filter_coconut_csv.py
```

Bước 2: Lọc dữ liệu từ tệp đầu ra tạo ở bước 1, mục tiêu tạo đầu ra là tệp .filtered.smi với 2 cột: SMILES và Name. với các chất đã qua các bộ lọc sau:
- Chất không có chứa các nguyên tố: 

```bash
python filter/filter_coconut_smi.py
```

Bước 3: Tách thành 3 tệp: train, valid, test
```bash
python filter/split_data.py data/pre-training/coconut/raw/coconut-complete-10-2024.filtered.smi --output /data/pre-training/coconut/
```

Bước 4: Preprocess

Bước 5: Train

Bước 6: Generate




