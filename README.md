# AutoFF

## Usage
```bash
git clone https://github.com/princesswinnie1122/AutoFF.git
```
Simply run "main.py" to call out the widget.

## Complete
- 可切割分子並在資料庫中搜尋片段 (過程 < 1分鐘)
- 生成腳本及計算所需檔案，視覺化確認結構無誤：[測試結果](https://docs.google.com/spreadsheets/d/1iIg0zPHXOcnBo-zzWsD-uY-R4ekoL3AcCz-gUzkuS94/edit?usp=sharing) (更新中)
- 判斷BB Not found程度，視情況要求使用者自行輸入、加入檔案
- 生成檔案前判斷分子資訊及檔案是否齊全，提醒缺少的部分
- 判斷SMILES是否前後一致，允許手動編輯腳本 (暫時解決假設接點都為Br的問題)


## To Do
- [ ] 簡化使用者新增BB拓樸流程 (用BRICS切割格式分類？)
- [ ] 解決特例 
  - EX1: Y17 的 N 會被切開 → 另外判斷 [5*]？
  - EX2: X5 不屬於BRICS結構，切割時會報錯 → 判斷後改成使用者自行定義BB？
- [ ] 後續計算部分的自動化

## More about the project
1. 辨認結構並產生初始檔案的工具
2. Non-bonded parameters 計算與整合的工具
3. Dihedral fitting 與檢查的工具
4. routine 分析 MD 表現的工具




