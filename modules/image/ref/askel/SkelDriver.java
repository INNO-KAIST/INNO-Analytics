// Fiji AnalyzeSkeleton_ 을 대화 상자 없이 실행해 결과를 줄 단위로 찍는다 (시험 정답용 — 앱에 들어가지 않음).
// 빌드: git clone https://github.com/fiji/AnalyzeSkeleton → javac -cp ij.jar -d out SkelDriver.java <AnalyzeSkeleton 의 6 파일>
// 실행: java -cp ij.jar:out SkelDriver <skeleton.tif> <prune 0..3> <pruneEnds 0/1> <lsp 0/1> <pw> <ph> <pd> [원본.tif]
import ij.*; import ij.measure.Calibration; import sc.fiji.analyzeSkeleton.*; import java.util.*;
public class SkelDriver {
  public static void main(String[] args) throws Exception {
    if (args[0].equals("-list")) { int k = 0; for (String line : java.nio.file.Files.readAllLines(java.nio.file.Paths.get(args[1]))) { if (line.trim().isEmpty()) continue; System.out.println("CASE " + (k++)); one(line.trim().split("\t")); } System.exit(0); }
    one(args); System.exit(0);
  }
  static void one(String[] a) {
    ImagePlus imp = IJ.openImage(a[0]);
    Calibration c = imp.getCalibration(); c.pixelWidth = Double.parseDouble(a[4]); c.pixelHeight = Double.parseDouble(a[5]); c.pixelDepth = Double.parseDouble(a[6]);
    ImagePlus orig = a.length > 7 ? IJ.openImage(a[7]) : null; if (orig != null) orig.setCalibration(c);
    AnalyzeSkeleton_ s = new AnalyzeSkeleton_(); s.setup("", imp);
    SkeletonResult r = s.run(Integer.parseInt(a[1]), a[2].equals("1"), a[3].equals("1"), orig, true, true);
    int n = r.getNumOfTrees(); System.out.println("TREES " + n);
    ArrayList<Double> sp = r.getShortestPathList(); double[][] sps = r.getSpStartPosition();
    for (int i = 0; i < n; i++) {
      StringBuilder b = new StringBuilder("T " + i + " " + r.getBranches()[i] + " " + r.getJunctions()[i] + " " + r.getEndPoints()[i] + " " + r.getJunctionVoxels()[i] + " " + r.getSlabs()[i]
        + " " + r.getAverageBranchLength()[i] + " " + r.getTriples()[i] + " " + r.getQuadruples()[i] + " " + r.getMaximumBranchLength()[i]);
      if (sp != null) b.append(" " + sp.get(i) + " " + sps[i][0] + " " + sps[i][1] + " " + sps[i][2]);
      System.out.println(b);
    }
    Graph[] g = r.getGraph();
    for (int i = 0; i < n; i++) for (Edge e : g[i].getEdges()) {
      Point p = e.getV1().getPoints().get(0), q = e.getV2().getPoints().get(0);
      System.out.println("E " + i + " " + e.getLength() + " " + p.x + " " + p.y + " " + p.z + " " + q.x + " " + q.y + " " + q.z + " " + e.getLength_ra() + " " + e.getSlabs().size());
    }
    ImageStack t = s.getResultImage(false); // 태그 영상: 끝 30, 줄기 127, 접점 70
    StringBuilder b = new StringBuilder("TAG");
    for (int z = 1; z <= t.getSize(); z++) { byte[] px = (byte[]) t.getPixels(z); for (int i = 0; i < px.length; i++) if (px[i] != 0) b.append(" " + ((z - 1) * px.length + i) + ":" + (px[i] & 255)); }
    System.out.println(b);
    if (a[3].equals("1")) { ImageStack l = s.getResultImage(true); StringBuilder c2 = new StringBuilder("LSPIMG"); // 가장 긴 최단 경로 영상: 값 96 인 화소
      if (l != null) for (int z = 1; z <= l.getSize(); z++) { byte[] px = (byte[]) l.getPixels(z); for (int i = 0; i < px.length; i++) if ((px[i] & 255) == 96) c2.append(" " + ((z - 1) * px.length + i)); }
      System.out.println(c2); }
  }
}
