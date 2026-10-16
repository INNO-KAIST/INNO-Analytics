/* Gwyddion 의 파일 가져오기 모듈로 원본을 읽어 GWY 로 쓴다 (시험 정답용).
   `gwyddion --convert-to-gwy` 는 GTK 를 초기화하지 않아 진행 창을 여는 모듈(ARDF 등)이 실패하므로 따로 만든다.
   빌드: gcc gwyconv.c -o gwyconv $(pkg-config --cflags --libs gwyddion)   실행: xvfb-run -a ./gwyconv 원본 출력.gwy */
#include <stdio.h>
#include <libgwyddion/gwyddion.h>
#include <libgwymodule/gwymodule.h>
#include <libprocess/gwyprocess.h>
#include <libgwydgets/gwydgets.h>
#include <app/gwyapp.h>
int main(int argc, char **argv) {
    GError *err = NULL; GwyContainer *c;
    gchar **dirs; gint i;
    if (argc != 3) { fprintf(stderr, "usage: gwyconv in out.gwy\n"); return 2; }
    gtk_init(&argc, &argv);
    gwy_widgets_type_init();
    dirs = gwy_app_settings_get_module_dirs();
    for (i = 0; dirs[i]; i++) if (g_getenv("GWYCONV_DEBUG")) fprintf(stderr, "module dir %s\n", dirs[i]);
    gwy_module_register_modules((const gchar**)dirs);
    c = gwy_file_load(argv[1], GWY_RUN_NONINTERACTIVE, &err);
    if (!c) { fprintf(stderr, "load failed: %s\n", err ? err->message : "?"); return 1; }
    if (!gwy_file_save(c, argv[2], GWY_RUN_NONINTERACTIVE, &err)) { fprintf(stderr, "save failed: %s\n", err ? err->message : "?"); return 1; }
    return 0;
}
