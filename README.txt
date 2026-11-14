iAM 全部计算源程序包
日期：2026-11-14

对应稿件：Further cases of the inductive Alperin–McKay condition
作者：Qiyue Tang；tangqy24@mails.tsinghua.edu.cn
对应范围：A、²A、B/C、整个幺幂 D 块、偶定义特征 D、F4。

本包包含本稿所用的全部 9 项有限计算源程序、实际输出、统一运行脚本、
逐项索引与文件校验。程序从已交付的最终源码包原样提取。
这是源程序包，不包含软件安装包、运行时或原始文献全文。
此前另行交付的 G2 与 ³D4 两族不属于这篇稿件。

一、程序与输出

checks/：
  check_a_hessian_exhaustive.py     Hessian 正规化子及 augmentation 边界
  check_root_models.py             quartic D4 与 E8 格模 2
  check_d4_actual_root_images.py   实际 D4 根自同构与中央化子
  o8_weyl_normalizers.g            实际 Weyl 群正规化子及局部普通字符
  o8_global_blocks.g               全局块、中央字符和外作用轨道
  quartic_independent.wl           quartic D4 的独立 Wolfram 计算

computations/：
  f4_independent_integer.py               F4 根与 Weyl 群的整数计算
  f4_small_field_source_round2.py         q=4,8 的实际混合环面源
  f4_small_field_parameters_round3.py     q=4,8 的实际有理半单参数

六个 Python 程序的实际结果为相邻的同名 JSON；两个 GAP 程序各有原始
stdout.txt 和 stderr.txt。Wolfram 实际输出为 checks/quartic_wolfram_result.json。
checks/o8_lattice_mod2.json 是格计算的额外记录；相同数据也由
check_root_models.py 的 o8_cover 部分计算。
calculation_index.json 列明每一项的输出、依赖及稿件章节。

二、完整性检查

解压后在本包目录运行：

  python verify_package.py

此命令检查文件 SHA-256，不执行数学计算。

三、运行 Python 与 GAP 计算

已使用的环境：Python 3.12.14，SymPy 1.14.0；GAP 4.14.0，CTblLib 1.3.9。
仅 check_root_models.py 和 f4_small_field_parameters_round3.py 需要 SymPy；
其余四个 Python 计算程序使用标准库。Python 依赖版本写在 requirements.txt。
GAP 脚本读取 CTblLib 自带的 2.O8+(2)、2.O8+(2).2 表及 fusion。
没有遗漏需要另行取得的自定义输入数据文件。

列出统一脚本包含的八项计算：

  python run_finite_checks.py --list

运行六个 Python 计算：

  python run_finite_checks.py --python-only

运行六个 Python 与两个 GAP 计算（gap 应在 PATH 中）：

  python run_finite_checks.py

指定 GAP 路径时：

  python run_finite_checks.py --gap /path/to/gap

若需指定 GAP 根目录，可另加 --gap-root /path/to/gap/root。
统一脚本将新日志写入 reproduced_outputs/；Python 计算还会更新相邻 JSON。
建议在解压副本中复跑，以保留交付时的原始输出。

四、独立 Wolfram 计算

第九项不在上述 Python/GAP 统一脚本中。请在新的 Wolfram kernel 中，
从本包目录执行：

  Get["checks/quartic_independent.wl"]

如已有 wolframscript，命令行可用：

  wolframscript -code 'Print[Get["checks/quartic_independent.wl"]]'

该源文件以 Association 作为返回值，命令行示例显式打印它。
记录的结果是 |W|=16、|R|=96，六个反射，以及六个阶 16 的阿贝尔中央化子。

五、证据范围

本次打包核对了全部源程序和输出的存在、JSON 结构、GAP 完成标记、
最终稿源码的一致性与文件完整性；未把打包检查记为一次新的计算复跑。
包内输出来自此前实际执行的计算。有限结构计算的使用范围由正文确定，
普通字符、块关系以及一般 iAM 论证仍见论文证明。
