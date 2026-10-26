# 新校正windowの独立検証

全210点のbest-speciesを、解析的な区間分割とは別の4核種closest-pullループで照合しました。同一pullのtie、幅差による分離した2区間、同じ幅の例も検証します。

新TPC窓と新best窓の信号・背景積分をraw dE/dx座標の独立Simpson積分と比較し、元gateの純度も再現します。全検証は校正を変更しません。

C++ runtime TSVが指定された場合は、実provider出力のcenter・width・校正適用flagをPythonと照合します。詳細・source SHAはvalidation.jsonに保存します。
