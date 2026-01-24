#include "turbo_encode.h"

#include <cstdlib>
#include <setjmp.h>

#include "jpeglib.h"

struct JpegError {
  jpeg_error_mgr pub;
  jmp_buf setjmp_buffer;
};

static void turbo_error_exit(j_common_ptr cinfo) {
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

bool encode_with_scans_turbo(const std::vector<uint8_t>& input_data,
                             const std::vector<ScanLine>& scans,
                             bool arithmetic,
                             std::vector<uint8_t>* out_bytes) {
  if (input_data.empty()) return false;

  jpeg_decompress_struct dinfo{};
  jpeg_compress_struct cinfo{};
  JpegError err{};
  jpeg_error_mgr* err_ptr = jpeg_std_error(&err.pub);
  err.pub = *err_ptr;
  err.pub.error_exit = turbo_error_exit;
  bool dinfo_created = false;
  bool cinfo_created = false;
  unsigned char* outbuffer = nullptr;
  if (setjmp(err.setjmp_buffer)) {
    if (cinfo_created) jpeg_destroy_compress(&cinfo);
    if (outbuffer) free(outbuffer);
    if (dinfo_created) jpeg_destroy_decompress(&dinfo);
    return false;
  }

  dinfo.err = &err.pub;
  jpeg_create_decompress(&dinfo);
  dinfo_created = true;
  jpeg_mem_src(&dinfo, const_cast<unsigned char*>(input_data.data()),
               static_cast<unsigned long>(input_data.size()));
  jpeg_read_header(&dinfo, TRUE);
  jvirt_barray_ptr* coef_arrays = jpeg_read_coefficients(&dinfo);

  cinfo.err = &err.pub;
  jpeg_create_compress(&cinfo);
  cinfo_created = true;

  unsigned long outsize = 0;
  jpeg_mem_dest(&cinfo, &outbuffer, &outsize);
  jpeg_copy_critical_parameters(&dinfo, &cinfo);
  cinfo.optimize_coding = !arithmetic;
  cinfo.arith_code = arithmetic;
  std::vector<jpeg_scan_info> scan_infos;
  if (!scans.empty()) {
    scan_infos.reserve(scans.size());
    for (const auto& scan : scans) scan_infos.push_back(to_jpeg_scan_info(scan));
    cinfo.scan_info = scan_infos.data();
    cinfo.num_scans = static_cast<int>(scan_infos.size());
  }

  jpeg_write_coefficients(&cinfo, coef_arrays);
  jpeg_finish_compress(&cinfo);

  if (out_bytes) out_bytes->assign(outbuffer, outbuffer + outsize);

  jpeg_destroy_compress(&cinfo);
  free(outbuffer);
  jpeg_finish_decompress(&dinfo);
  jpeg_destroy_decompress(&dinfo);
  return true;
}
