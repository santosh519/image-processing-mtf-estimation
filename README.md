# Modulation Transfer Function (MTF) Analysis of IKONOS Panchromatic Imagery

This repository contains a complete MATLAB workflow for calculating the **Modulation Transfer Function (MTF)** from a high-resolution satellite image. The analysis is based on a panchromatic band of **IKONOS imagery acquired on August 1, 2005**. The workflow includes edge detection, Edge Spread Function (ESF) oversampling and smoothing, Line Spread Function (LSF) generation, and final MTF derivation.

> **Note**: The image data used in this project is licensed and cannot be redistributed. However, the panchromatic IKONOS data used in this project can be downloaded from USGS Earth Explorer, European Space Agency (ESA) Archive or Maxar Technologies for authorized users.

---

## 📁 Project Structure

- `MTF_Project_IKONOS.m`: The main MATLAB script that executes the full MTF analysis pipeline.
- `MTF Project Report.pdf`: Detailed report explaining the methodology, results, and interpretation of the MTF analysis.
- `README.md`: This documentation file.

---

## 🔍 Analysis Workflow

1. **Image Loading and Preprocessing**
   - Reads a high-resolution TIFF image using `readgeoraster`.
   - Visualizes the image and applies binary masking to remove zero-value pixels.

2. **Region of Interest (ROI) Selection**
   - A rectangular ROI is selected using map coordinates and converted to pixel coordinates.
   - The image is cropped to include the ROI edge.

3. **Edge Spread Function (ESF) Extraction**
   - Fits a cubic spline to each scanline in the ROI.
   - Determines sub-pixel inflection points to estimate the edge location.
   - Projects pixel coordinates perpendicular to the edge direction to form a raw ESF.

4. **Savitzky-Golay Smoothing**
   - Applies a Savitzky-Golay filter with a sliding window and polynomial fit to generate a clean, uniformly sampled ESF.

5. **Signal-to-Noise Ratio (SNR) Estimation**
   - Calculates mean intensities in background and foreground regions.
   - Computes SNR and overlays results on ESF plot.

6. **Line Spread Function (LSF) and Full Width at Half Maximum (FWHM)**
   - Derives LSF by taking the derivative of the ESF.
   - Smooths LSF and computes the FWHM using interpolation.

7. **Modulation Transfer Function (MTF) Computation**
   - Applies FFT on the LSF.
   - Computes spatial frequency axis and plots normalized MTF against cycles per pixel.

---

## 📊 Output Visualizations

- ROI Cropped Edge Profile
- ESF (raw and smoothed)
- LSF with FWHM
- MTF Curve with Nyquist Frequency Annotation

---

## ⚙️ Requirements

- MATLAB R2021b or later
- Image Processing Toolbox
- Mapping Toolbox
- Curve Fitting Toolbox (for `fit()`)

---

## 🔒 Data Disclaimer

This code was developed and tested using **IKONOS panchromatic imagery (1m resolution)**, acquired on **August 1, 2005**. Due to licensing restrictions, the image cannot be shared publicly. Interested users should obtain access through authorized vendors.

## 👤 Author & Contact
Santosh Adhikari

Email: santadh2015@gmail.com

GitHub: @santosh519

Thank you for reviewing! Feedback and contributions are welcome.

Once obtained, place the `.tif` file in the appropriate folder as specified in the script under:

```matlab
ParentFolderPath = '...';
