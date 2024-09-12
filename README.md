# AutoFF

## Usage
'''=bash
git clone https://github.com/princesswinnie1122/AutoFF.git
'''
Simply run the "main.py" program to call out the widget.

## Complete
- 可切割分子並在資料庫中搜尋片段 (過程 < 1分鐘)
- 生成計算所需參數檔，視覺化確認結構無誤：測試結果 (更新中)
- 生成檔案後會比對 SMILES，若不同則提醒使用者自行編輯腳本

## To Do
[] 簡化使用者新增BB拓樸流程 (用BRICS切割格式分類？)
[] 解決特例 
    - EX1: Y17 的 N 會被切開 → 另外判斷 [5*]？
    - EX2: X5 不屬於BRICS結構，切割時會報錯 → 判斷後改成使用者自行定義BB？
[] 後續計算部分的自動化
[] 目前假設所有接點都是 Br，無法適應接點不同造成同BB不同結構的情形 

## More about the project
1. 辨認結構並產生初始檔案的工具
2. Non-bonded parameters 計算與整合的工具
3. Dihedral fitting 與檢查的工具
4. routine 分析 MD 表現的工具