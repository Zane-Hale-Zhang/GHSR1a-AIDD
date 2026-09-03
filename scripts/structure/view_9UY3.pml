# ==========================================================
# GHSR1a 结构查看脚本 · 9UY3 (anamorelin 复合物, 2.52 Å)
#
# 用法（三选一）:
#   1. PyMOL 命令行执行:  @~/aidd/structures/view_9UY3.pml
#   2. PyMOL 菜单: File → Run Script → 选这个文件
#   3. 终端: pymol ~/aidd/structures/9UY3.cif -r ~/aidd/structures/view_9UY3.pml
# ==========================================================

# ---------- 加载 ----------
load ~/aidd/structures/9UY3.cif, ghsr
hide everything, ghsr

# ---------- 受体本体（链 R），其余 G 蛋白/纳米抗体不显示 ----------
show cartoon, ghsr and chain R
set cartoon_fancy_helices, 1
color purple, ghsr and chain R

# ---------- 共晶配体 anamorelin (A1ESD) ----------
select ligand, ghsr and resn A1ESD
show sticks, ligand
color red, ligand
set stick_radius, 0.28, ligand

# ---------- 口袋残基（配体 5 Å 内）----------
select pocket, ghsr and chain R and byres (ligand around 5)
show sticks, pocket
color gray50, pocket
set stick_transparency, 0.25, pocket

# ---------- 只显示口袋残基的编号，看清是哪几个 ----------
set label_size, 14
set label_color, black
label pocket and name CA, "%s%s" % (resn, resi)

# ---------- 视角与渲染 ----------
orient ghsr and chain R
zoom ligand, 14

bg_color white
set antialias, 2
set ray_shadows, 0
set cartoon_transparency, 0.2
set sphere_quality, 2

# ---------- 出图 ----------
ray 1400, 1000
png ~/aidd/structures/ghsr1a_overview.png

print ""
print "  ─────────────────────────────────────────"
print "  紫色卡通 = GHSR1a 受体 (链 R)"
print "  红色棍状 = anamorelin (共晶配体)"
print "  灰色棍状 = 口袋残基（配体 5 Å 内）"
print "  图已保存: ~/aidd/structures/ghsr1a_overview.png"
print "  ─────────────────────────────────────────"
print ""
print "  试着在命令行里敲这些，看会发生什么:"
print "    hide cartoon                 # 只留配体和口袋"
print "    show surface, chain R        # 看口袋的形状"
print "    set surface_transparency,0.5 # 表面半透明，能看见里面的配体"
print "    select pocket and name CA    # 选中口袋残基，看左下角列出了哪些"
print "    zoom ligand, 8               # 放大看结合模式"
print ""
