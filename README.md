# Structural-Sensitivity-Analysis

ReactionNetwork.py
reactionのリストを与えてネットワークを作成．
限局構造の探索，グラフ指数の計算を行う

find_limitset
時間がかかるので非推奨．A行列を乱数生成してその度限局集合を計算し，現れたすべての限局集合を返す関数．
S行列の0成分を判定する閾値は1.0e-10.

find_limitset_meansmat
S行列を複数回生成してその平均をとる．
0成分を判定する閾値は1.0e-10だが，閾値周辺の値がある (1.0e-8~1.0e-10) と警告を返す．
大きいネットワークだと警告が出る
