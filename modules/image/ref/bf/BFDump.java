// Bio-Formats 로 파일을 읽어 series 별 크기·형식·보정·채널 이름을 찍고, 화소를 float32 LE 로 쏟는다 (시험 정답용 — 앱에 들어가지 않음)
// 빌드: javac -cp bioformats_package.jar -d out BFDump.java   실행: java -cp bioformats_package.jar:out BFDump <파일> [<출력 폴더>] [series 목록 "0,3"]
import loci.formats.*; import loci.formats.meta.*; import loci.common.services.*; import loci.formats.services.*; import ome.units.UNITS;
import java.io.*; import java.nio.*; import java.util.*;
public class BFDump {
  public static void main(String[] a) throws Exception {
    loci.common.DebugTools.setRootLevel("ERROR");
    ServiceFactory f = new ServiceFactory(); OMEXMLService svc = f.getInstance(OMEXMLService.class); IMetadata md = svc.createOMEXMLMetadata();
    ImageReader r = new ImageReader(); r.setMetadataStore(md); r.setFlattenedResolutions(false); r.setId(a[0]);
    String out = a.length > 1 ? a[1] : null; Set<Integer> want = null;
    if (a.length > 2) { want = new HashSet<>(); for (String s : a[2].split(",")) want.add(Integer.parseInt(s)); }
    System.out.println("FORMAT " + r.getFormat() + " SERIES " + r.getSeriesCount());
    for (int s = 0; s < r.getSeriesCount(); s++) {
      r.setSeries(s);
      String px = md.getPixelsPhysicalSizeX(s) == null ? "-" : "" + md.getPixelsPhysicalSizeX(s).value(UNITS.MICROMETER);
      String py = md.getPixelsPhysicalSizeY(s) == null ? "-" : "" + md.getPixelsPhysicalSizeY(s).value(UNITS.MICROMETER);
      String pz = md.getPixelsPhysicalSizeZ(s) == null ? "-" : "" + md.getPixelsPhysicalSizeZ(s).value(UNITS.MICROMETER);
      StringBuilder ch = new StringBuilder(); for (int c = 0; c < md.getChannelCount(s); c++) ch.append(c == 0 ? "" : "|").append(md.getChannelName(s, c));
      System.out.println("S " + s + " X " + r.getSizeX() + " Y " + r.getSizeY() + " Z " + r.getSizeZ() + " C " + r.getSizeC() + " T " + r.getSizeT() + " eff " + r.getEffectiveSizeC() + " rgb " + r.getRGBChannelCount()
        + " type " + FormatTools.getPixelTypeString(r.getPixelType()) + " bits " + r.getBitsPerPixel() + " order " + r.getDimensionOrder() + " interleaved " + r.isInterleaved() + " little " + r.isLittleEndian()
        + " px " + px + " py " + py + " pz " + pz + " name " + md.getImageName(s).replace(' ', '_') + " ch " + ch.toString().replace(' ', '_'));
      if (out == null || (want != null && !want.contains(s))) continue;
      // 화소: 앱 순서 (c, z, t) → c + C·(z + Z·t), RGB 는 표본별 장
      int X = r.getSizeX(), Y = r.getSizeY(), rgb = r.getRGBChannelCount(), bpp = FormatTools.getBytesPerPixel(r.getPixelType()); int pt = r.getPixelType();
      try (DataOutputStream o = new DataOutputStream(new BufferedOutputStream(new FileOutputStream(out + "/s" + s + ".f32")))) {
        for (int t = 0; t < r.getSizeT(); t++) for (int z = 0; z < r.getSizeZ(); z++) for (int c = 0; c < r.getEffectiveSizeC(); c++) {
          byte[] b = r.openBytes(r.getIndex(z, c, t)); ByteBuffer bb = ByteBuffer.wrap(b).order(r.isLittleEndian() ? ByteOrder.LITTLE_ENDIAN : ByteOrder.BIG_ENDIAN);
          for (int k = 0; k < rgb; k++) for (int i = 0; i < X * Y; i++) {
            int idx = r.isInterleaved() ? (i * rgb + k) : (k * X * Y + i); float v;
            switch (pt) { case FormatTools.UINT8: v = b[idx] & 255; break; case FormatTools.INT8: v = b[idx]; break; case FormatTools.UINT16: v = bb.getShort(idx * 2) & 0xffff; break; case FormatTools.INT16: v = bb.getShort(idx * 2); break;
              case FormatTools.UINT32: v = (float)(bb.getInt(idx * 4) & 0xffffffffL); break; case FormatTools.INT32: v = bb.getInt(idx * 4); break; case FormatTools.FLOAT: v = bb.getFloat(idx * 4); break; case FormatTools.DOUBLE: v = (float) bb.getDouble(idx * 8); break; default: throw new RuntimeException("type " + pt); }
            o.writeInt(Integer.reverseBytes(Float.floatToIntBits(v))); }
        }
      }
    }
    r.close();
  }
}
