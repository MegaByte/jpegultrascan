#include "mozjpeg_encode.h"

#include <cstdlib>
#include <setjmp.h>

#include "jpeglib.h"

struct JpegError {
  jpeg_error_mgr pub;
  jmp_buf setjmp_buffer;
};

static void mozjpeg_error_exit(j_common_ptr cinfo) {
  JpegError* err = reinterpret_cast<JpegError*>(cinfo->err);
  (*cinfo->err->output_message)(cinfo);
  longjmp(err->setjmp_buffer, 1);
}

static jpeg_scan_info to_jpeg_scan_info(const ScanLine& scan) {
  jpeg_scan_info info{};
  info.comps_in_scan = static_cast<int>(scan.components.size());
  for (size_t i = 0; i < scan.components.size(); ++i) {
    info.component_index[i] = scan.components[i];
  }
  info.Ss = scan.Ss;
  info.Se = scan.Se;
  info.Ah = scan.Ah;
  info.Al = scan.Al;
  return info;
}

static bool encode_with_scans_mozjpeg_internal(j_decompress_ptr dinfo,
                                               jvirt_barray_ptr* coef_arrays,
                                               const std::vector<ScanLine>& scans,
                                               bool minimal_scan,
                                               std::vector<uint8_t>* out_bytes) {
  jpeg_compress_struct cinfo{};
  JpegError cerr{};
  cinfo.err = jpeg_std_error(&cerr.pub);
  cerr.pub.error_exit = mozjpeg_error_exit;
  if (setjmp(cerr.setjmp_buffer)) {
    jpeg_destroy_compress(&cinfo);
    return false;
  }
  jpeg_create_compress(&cinfo);

  unsigned char* outbuffer = nullptr;
  unsigned long outsize = 0;
  jpeg_mem_dest(&cinfo, &outbuffer, &outsize);
  jpeg_copy_critical_parameters(dinfo, &cinfo);
  cinfo.optimize_coding = FALSE;
  cinfo.arith_code = FALSE;

  std::vector<jpeg_scan_info> scan_infos;
  if (!minimal_scan && !scans.empty()) {
    scan_infos.reserve(scans.size());
    for (const auto& scan : scans) scan_infos.push_back(to_jpeg_scan_info(scan));
    cinfo.scan_info = scan_infos.data();
    cinfo.num_scans = static_cast<int>(scan_infos.size());
    cinfo.progressive_mode = TRUE;
    jpeg_c_set_bool_param(&cinfo, JBOOLEAN_OPTIMIZE_SCANS, FALSE);
  } else if (minimal_scan) {
    cinfo.scan_info = nullptr;
    cinfo.num_scans = 0;
    cinfo.progressive_mode = FALSE;
  }

  jpeg_write_coefficients(&cinfo, coef_arrays);
  jpeg_finish_compress(&cinfo);

  if (out_bytes) out_bytes->assign(outbuffer, outbuffer + outsize);

  jpeg_destroy_compress(&cinfo);
  free(outbuffer);
  return true;
}

bool encode_with_scans_mozjpeg_from_coeffs(j_decompress_ptr dinfo,
                                           jvirt_barray_ptr* coef_arrays,
                                           const std::vector<ScanLine>& scans,
                                           bool /*arithmetic*/,
                                           bool minimal_scan,
                                           std::vector<uint8_t>* out_bytes) {
  if (!dinfo || !coef_arrays) return false;
  return encode_with_scans_mozjpeg_internal(dinfo, coef_arrays, scans,
                                            minimal_scan, out_bytes);
}

bool encode_with_scans_mozjpeg(const std::vector<uint8_t>& input_data,
                               const std::vector<ScanLine>& scans,
                               bool /*arithmetic*/,
                               std::vector<uint8_t>* out_bytes) {
  if (input_data.empty()) return false;

  jpeg_decompress_struct dinfo{};
  JpegError derr{};
  dinfo.err = jpeg_std_error(&derr.pub);
  derr.pub.error_exit = mozjpeg_error_exit;
  if (setjmp(derr.setjmp_buffer)) {
    jpeg_destroy_decompress(&dinfo);
    return false;
  }
  jpeg_create_decompress(&dinfo);
  jpeg_mem_src(&dinfo, const_cast<unsigned char*>(input_data.data()),
               static_cast<unsigned long>(input_data.size()));
  jpeg_read_header(&dinfo, TRUE);
  jvirt_barray_ptr* coef_arrays = jpeg_read_coefficients(&dinfo);

  bool ok = encode_with_scans_mozjpeg_internal(&dinfo, coef_arrays, scans, false, out_bytes);
  jpeg_finish_decompress(&dinfo);
  jpeg_destroy_decompress(&dinfo);
  return ok;
}
