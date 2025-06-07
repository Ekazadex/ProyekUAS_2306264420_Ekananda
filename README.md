# Implementasi Metode Runge-Kutta Orde 4 untuk Penyelesaian Persamaan Diferensial Biasa

- **Nama**: Ekananda Zhafif Dean
- **NPM**: 2306264420

## Deskripsi Program

Program ini merupakan implementasi komprehensif dari **Metode Runge-Kutta Orde 4** untuk menyelesaikan persamaan diferensial biasa (ODE) dengan berbagai fitur analisis dan perbandingan. Program dikembangkan dalam bahasa C++ dengan pendekatan berorientasi objek untuk memberikan solusi yang akurat dan efisien.

### Fitur Utama

#### 1. **Metode Numerik yang Diimplementasikan**
- **Runge-Kutta Orde 4 (RK4)** - Metode utama dengan akurasi tinggi
- **Metode Euler** - Metode dasar untuk perbandingan
- **Runge-Kutta Orde 2 (Heun)** - Metode dengan akurasi sedang
- **Runge-Kutta Orde 3** - Metode dengan akurasi tinggi alternatif

#### 2. **Analisis Error dan Akurasi**
- Perhitungan **error absolut** dan **error relatif**
- Perbandingan dengan solusi analitik
- Analisis maksimum error
- Validasi akurasi metode

#### 3. **Analisis Konvergensi**
- Uji konvergensi dengan berbagai ukuran step
- Perhitungan orde konvergensi
- Analisis stabilitas numerik
- Perbandingan waktu komputasi

#### 4. **Kasus Uji Komprehensif**
- **Kasus 1**: Pertumbuhan Eksponensial (`dy/dx = y`)
- **Kasus 2**: Peluruhan Gaussian (`dy/dx = -2xy`)
- **Kasus 3**: Persamaan Linear Non-homogen (`dy/dx = x + y`)

#### 5. **Analisis Stabilitas**
- Uji stabilitas dengan berbagai parameter
- Analisis batas stabilitas metode
- Testing dengan persamaan kaku (stiff equations)

#### 6. **Output dan Visualisasi Data**
- Format output yang terstruktur dan mudah dibaca
- Export hasil ke file CSV
- Laporan analisis yang komprehensif
- Perbandingan performa antar metode

### Keunggulan Program

1. **Akurasi Tinggi**: Metode RK4 memberikan akurasi orde 4 dengan error yang sangat kecil
2. **Efisiensi Komputasi**: Optimasi kode untuk kecepatan eksekusi
3. **Fleksibilitas**: Dapat menangani berbagai jenis persamaan diferensial
4. **Analisis Mendalam**: Menyediakan berbagai metrik untuk evaluasi performa
5. **Portabilitas**: Kode C++ yang dapat dikompilasi di berbagai platform
