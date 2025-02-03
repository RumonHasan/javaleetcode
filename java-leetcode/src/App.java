import java.util.*;

public class App {
    public static void main(String[] args) throws Exception {
        // longestConsequtiveNums(new int[] { 100, 4, 200, 1, 3, 2 });
        // lexicologicalSmallestSwaps(new int[] { 1, 7, 6, 18, 2, 1 }, 3);
        // longestCommonPrefix(new String[] { "reflower", "flow", "flight" });
        // int[][] grid = {
        // { 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0 },
        // { 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 0, 0, 0 },
        // { 0, 1, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0 },
        // { 0, 1, 0, 0, 1, 1, 0, 0, 1, 0, 1, 0, 0 },
        // { 0, 1, 0, 0, 1, 1, 0, 0, 1, 1, 1, 0, 0 },
        // { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0 },
        // { 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 0, 0, 0 },
        // { 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0 }
        // };
        // maxAreaOfIsland(grid);
        // mostCompetitiveArray(new int[] { 3, 5, 2, 6 }, 2);
//        // smallestSubstring(new String("cbacdcbc"));
//        coinChange(new int[] { 1, 2, 5 }, 11);
//        int[][] grid = {
//                { 0, 2, 1, 0 },
//                { 4, 0, 0, 3 },
//                { 1, 0, 0, 4 },
//                { 0, 3, 2, 0 }
//        };

        // findMaxFish(grid);
        // pushDominoes(new String(".L.R...LR..L.."));
        // optimalStringPartition(new String("abacaba"));
        //maximumSubarraySum(new int[] { 1, 5, 4, 2, 9, 9, 9 }, 3);
//        maximumChunksSubarraySorted(new int[] {1,0,2,3,4});
//        longestStringChain(new String[] {"a","b","ba","bca","bda","bdca"});
//        int[][] matrix = {
//                {0, 0, 0, 0},
//                {1, 0, 1, 0},
//                {0, 1, 1, 0},
//                {0, 0, 0, 0}
//        };
//
//        numberOfEnclaves(matrix);
//        int[][] farmland = {
//                {1, 0, 0},
//                {0, 1, 1},
//                {0, 1, 1}
//        };
//        groupFarmlands(farmland);
        longestMonotonicSubarray(new int [] {1,4,3,3,2});
        char[][] board = {
                {'X', 'X', 'X', 'X'},
                {'X', 'O', 'O', 'X'},
                {'X', 'X', 'O', 'X'},
                {'X', 'O', 'X', 'X'}
        };
        solve(board);


    }

    // longest consequtive sum
    public static int longestConsequtiveNums(int[] nums) {
        HashSet<Integer> numSet = new HashSet<Integer>();
        long maxLen = 0;
        for (int num : nums) {
            numSet.add(num);
        }

        // checking for consequtives
        for (int j = 0; j < nums.length; j++) {
            int curr = nums[j];
            if (!numSet.contains(curr + 1)) {
                long counter = 1;
                long prevNum = (long) curr - 1; // type casting to match the int number
                while (numSet.contains((int) prevNum)) {
                    prevNum--;
                    counter++;
                }
                maxLen = Math.max(maxLen, counter);
            }
        }
        return (int) maxLen;
    }

    // dividing into groups then swapping around the values
    public static int[] lexicologicalSmallestSwaps(int[] nums, int limit) {

        Integer[] sortedNums = Arrays.stream(nums).boxed().toArray(Integer[]::new);
        Arrays.sort(sortedNums, (a, b) -> a - b);

        int[] result = new int[nums.length];
        List<List<Integer>> groups = new ArrayList<>();
        List<Integer> currentGroup = new ArrayList<>();
        Map<Integer, Integer> groupMap = new HashMap<>();
        int groupIndex = 0;

        // Group formation
        for (int i = 0; i < sortedNums.length; i++) {
            int currNum = sortedNums[i];

            if (currentGroup.isEmpty() || // if its empty or has element then add
                    Math.abs(currentGroup.get(currentGroup.size() - 1) - currNum) <= limit) {
                currentGroup.add(currNum);
                groupMap.put(currNum, groupIndex);
            } else {
                groups.add(new ArrayList<>(currentGroup));
                currentGroup.clear();
                currentGroup.add(currNum);
                groupIndex++;
                groupMap.put(currNum, groupIndex);
            }

            // Add last group
            if (i == sortedNums.length - 1) {
                groups.add(new ArrayList<>(currentGroup));
            }
        }
        // Populate result
        for (int i = 0; i < nums.length; i++) {
            int currNum = nums[i];
            int groupIdx = groupMap.get(currNum);
            result[i] = groups.get(groupIdx).remove(0);
        }

        return result;
    }

    // basic common longest prefix
    public static String longestCommonPrefix(String[] strs) {
        StringBuilder finalPref = new StringBuilder("");
        StringBuilder prefix = new StringBuilder("");
        char[] checkWord = strs[0].toCharArray();
        String[] sliceStrs = Arrays.copyOfRange(strs, 1, strs.length);
        // checking prefix against each word
        for (char checkWordChar : checkWord) {
            prefix.append(checkWordChar);
            int counter = 0;
            for (int i = 0; i < sliceStrs.length; i++) {
                String currWord = sliceStrs[i];
                if (currWord.startsWith(prefix.toString())) {
                    counter++;
                } else {
                    break;
                }
            }
            if (counter == sliceStrs.length && finalPref.length() < prefix.length()) {
                finalPref = new StringBuilder(prefix);
            } else {
                break;
            }

        }
        return finalPref.toString();
    }

    // question for max island

    public static int dfs(int row, int col, int[][] grid, int mainRow, int mainCol) {
        if (row < 0 || col < 0 || row >= mainRow || col >= mainCol || grid[row][col] == 0) {
            return 0;
        }
        ;

        grid[row][col] = 0;// updating the existing values to 0 to prevent recount
        // recursive checks for up down left and right
        return (1 + dfs(row + 1, col, grid, mainRow, mainCol) + dfs(row, col + 1, grid, mainRow, mainCol)
                + dfs(row - 1, col, grid, mainRow, mainCol) + dfs(row, col - 1, grid, mainRow, mainCol));
    }

    public static int maxAreaOfIsland(int[][] grid) {
        int maxArea = 0;
        int ROW = grid.length;
        int COL = grid[0].length;

        // from every starting cell running dfs to check
        for (int i = 0; i < ROW; i++) {
            for (int j = 0; j < COL; j++) {
                maxArea = Math.max(maxArea, dfs(i, j, grid, ROW, COL));
            }
        }
        return maxArea;
    }

    // getting the most competitive array
    public static int[] mostCompetitiveArray(int[] nums, int k) {

        List<Integer> stack = new ArrayList<>();
        int LEN = nums.length;

        for (int index = 0; index < nums.length; index++) {
            int currNum = nums[index];
            while (!stack.isEmpty() && stack.get(stack.size() - 1) > currNum && stack.size() + (LEN - index - 1) >= k) {
                stack.remove(stack.size() - 1); // checking how many u can pop
            }
            if (stack.size() < k) {
                stack.add(currNum);
            }
        }
        int[] resultStack = new int[stack.size()];
        Arrays.fill(resultStack, 0);
        // just for printing
        for (int i = 0; i < resultStack.length; i++) {
            resultStack[i] = stack.get(i);
        }
        return resultStack;
    }

    // getting the smallest subsequence
    public static String smallestSubstring(String s) {
        // variables to populate the stack and substring
        StringBuilder resultStr = new StringBuilder();
        List<Character> stack = new ArrayList<>();
        HashMap<Character, Integer> map = new HashMap<>();
        HashSet<Character> set = new HashSet<>();
        char[] str = s.toCharArray();

        // getting occurence
        for (char localChar : str) {
            if (map.containsKey(localChar)) {
                map.put(localChar, map.get(localChar) + 1);
            } else {
                map.put(localChar, 1);
            }
        }
        ;

        // populating array list
        for (int i = 0; i < str.length; i++) {
            char currChar = str[i];
            // reducing from the map
            map.put(currChar, map.get(currChar) - 1);
            if (!set.contains(currChar)) {
                // removing if the char is bigger and there is an option of removing and adding
                // later
                while (!stack.isEmpty() && stack.get(stack.size() - 1) > currChar
                        && map.get(stack.get(stack.size() - 1)) > 0) {
                    set.remove(stack.get(stack.size() - 1));
                    stack.remove(stack.size() - 1);
                }
                // adding the char if there is no option
                stack.add(currChar);
                set.add(currChar);
            }

        }
        for (int i = 0; i < stack.size(); i++) {
            resultStr.append(stack.get(i));
        }

        return resultStr.toString();
    }
    // "cbacdcbc"

    public static boolean canPlaceFlowers(int[] flowerbed, int n) {
        int counter = n;

        for (int index = 0; index < flowerbed.length; index++) {
            int curr = flowerbed[index];

            if (curr == 0) {
                if (index == 0 && flowerbed[index + 1] != 1) {
                    flowerbed[index] = 1;
                    counter--;
                }
                ;
                if (index == flowerbed.length - 1 && flowerbed[index - 1] != 1) {
                    flowerbed[index] = 1;
                    counter--;
                }
                ;
                if (index > 0 && index < flowerbed.length - 1 && flowerbed[index - 1] != 1
                        && flowerbed[index + 1] != 1) {
                    flowerbed[index] = 1;
                    counter--;
                }
            }
            ;

            if (counter == 0) {
                break;
            }
        }
        boolean check = counter == 0;
        return check;
    }

    // using dynamic programmming approach
    public static int coinChange(int[] coins, int amount) {

        int[] dpArray = new int[amount + 1];
        Arrays.fill(dpArray, amount + 1); // integer max to predefine array
        dpArray[0] = 0;// 0 coins for 0 amount

        // dynamic programming iteration
        for (int amountIndex = 1; amountIndex <= amount; amountIndex++) {
            int currAmount = amountIndex;
            for (int i = 0; i < coins.length; i++) {
                int currCoin = coins[i];

                // check for dp iteration if only coin is smaller
                if (currCoin <= currAmount) {
                    int checkAmountDiff = currAmount - currCoin;
                    // checking with the existing dpArray value or adding with new coin difference
                    dpArray[amountIndex] = Math.min(dpArray[amountIndex], dpArray[checkAmountDiff] + 1);
                }
            }
        }
        return dpArray[amount] > amount ? -1 : dpArray[amount];
    };

    // fish in a pond dfs problem

    // dfs call to check for fish grid by marking the visited cell as 0 when visited
    public static int findMaxFishDfs(int row, int col, int[][] grid, int ROW, int COL) {
        // edge case
        if (row < 0 || col < 0 || row >= ROW || col >= COL || grid[row][col] == 0) {
            return 0;
        }
        int localMaxFishCount = 0;
        localMaxFishCount += grid[row][col]; // includes the current count of fish
        grid[row][col] = 0;
        // recursive call for the adjacent cells
        int left = findMaxFishDfs(row - 1, col, grid, ROW, COL);
        int right = findMaxFishDfs(row + 1, col, grid, ROW, COL);
        int top = findMaxFishDfs(row, col - 1, grid, ROW, COL);
        int bottom = findMaxFishDfs(row, col + 1, grid, ROW, COL);

        return (left + right + top + bottom + localMaxFishCount);
    }

    public static int findMaxFish(int[][] grid) {
        int maxFishCount = 0;
        int ROW = grid.length;
        int COL = grid[0].length;

        for (int i = 0; i < ROW; i++) {
            for (int j = 0; j < COL; j++) {
                int currGridVal = grid[i][j];
                if (currGridVal > 0) {
                    maxFishCount = Math.max(maxFishCount, findMaxFishDfs(i, j, grid, ROW, COL));
                }

            }
        }

        return maxFishCount;
    }

    // pushing dominoes problem
    public static String pushDominoes(String dominoes) {
        StringBuilder result = new StringBuilder();
        char[] domChars = dominoes.toCharArray();
        char[] resultChars = new char[domChars.length];
        Arrays.fill(resultChars, ' ');
        // for storing the forces of motion for left and right
        int[] left = new int[domChars.length];
        Arrays.fill(left, 0);
        int[] right = new int[domChars.length];
        Arrays.fill(right, 0);

        int charForce = 0;

        // right
        for (int i = 0; i < domChars.length; i++) {
            char currChar = domChars[i];
            if (currChar == 'R') {
                charForce = domChars.length;
            } else if (currChar == '.') {
                charForce--;
            } else {
                charForce = 0;
            }
            right[i] = Math.max(0, charForce);
        }
        ;
        // left
        for (int i = domChars.length - 1; i >= 0; i--) {
            char currChar = domChars[i];
            if (currChar == 'L') {
                charForce = domChars.length;
            } else if (currChar == '.') {
                charForce--;
            } else {
                charForce = 0;
            }
            left[i] = Math.max(0, charForce);
        }
        ;
        // populating the array
        for (int i = 0; i < resultChars.length; i++) {
            int leftVal = left[i];
            int rightVal = right[i];

            if (leftVal > rightVal) {
                resultChars[i] = 'L';
            } else if (rightVal > leftVal) {
                resultChars[i] = 'R';
            } else {
                resultChars[i] = '.';
            }

        }
        ;
        for (int i = 0; i < resultChars.length; i++) {
            result.append(resultChars[i]);
        }
        return result.toString();
    }

    // ".L.R...LR..L.."

    public static int optimalStringPartition(String s) {
        char[] str = s.toCharArray();
        HashSet<Character> set = new HashSet<>();
        int partitionCounter = 0;

        for (int i = 0; i < str.length; i++) {
            char currChar = str[i];

            // if set has a char already then increase counter and clear set
            if (set.contains(currChar)) {
                partitionCounter++;
                set.clear();
            }

            // add to set
            set.add(currChar); // for next partition
            if (i == str.length - 1) {
                partitionCounter++;
            }

        }

        return partitionCounter;
    }

    // checking for max sum when subarray length is distinct
    public static long maximumSubarraySum(int[] nums, int k) {
        long maxSum = 0;
        Map<Integer, Integer> map = new HashMap<>(); // Using Integer is sufficient for counts
        long localSum = 0;

        // Setting initial subarray count
        for (int i = 0; i < k; i++) {
            localSum += (long) nums[i]; // Cast to long before adding
            map.merge(nums[i], 1, Integer::sum);
        }

        maxSum = map.size() == k ? localSum : 0; // Initial max Sum
        int start = 0;

        for (int i = k; i < nums.length; i++) {
            int currNum = nums[i];
            int startNum = nums[start];

            localSum += currNum;
            localSum -= startNum;

            // here in case of map in java.. the key entry is removed if it returns null
            // since it returns null
            map.merge(startNum, -1, (oldValue, one) -> {
                int newValue = oldValue + one;
                return newValue == 0 ? null : newValue;
            });

            map.merge(currNum, 1, Integer::sum);

            if (map.size() == k) {
                maxSum = Math.max(localSum, maxSum);
            }
            start++;
        }

        return maxSum;

    }

    // question to find the maximum chunks that can be reconnected to make a new subarray
    public static int maximumChunksSubarraySorted(int[] arr){
        int maxChunks = 1; // by default there is always a single partition
        List<Integer> prefixMax = new ArrayList<>();
        int[] suffixMin = new int[arr.length];
        Arrays.fill(suffixMin, 0);
        suffixMin[suffixMin.length - 1] = arr[arr.length - 1]; // initialising with the last element for min calculation

        // getting the prefix max from left to right
        for(int i = 0; i < arr.length; i++){
            if(i == 0){
                prefixMax.add(arr[i]);
            }else{
                if(prefixMax.getLast() < arr[i]){
                    prefixMax.add(arr[i]);
                }else{
                    prefixMax.add(prefixMax.getLast());
                }
            }
        }
        // getting the suffix min so will start from the right side biggest to smallest left side
        for(int i = arr.length - 2; i >= 0; i--){
            int currNum = arr[i];
            suffixMin[i] = Math.min(currNum, suffixMin[i + 1]);
        };
        // getting the max chunks
        for(int i = 1; i < arr.length; i++){
            int rightSmallest = suffixMin[i];
            int currElement = prefixMax.get(i - 1);
            // if the right smallest is bigger than left side then its a possible partition
            if(currElement <= rightSmallest){
                maxChunks++;
            }
        }

        return maxChunks;
    }


    // getting longest string chain
    // note: dfs should be applied on every word as the starting point of the chain
    public static int longestStringChain(String[] words){
        int maxChain = 0;
        HashMap<String, Integer> map = new HashMap<>();
        String[] sortedWords = Arrays.stream(words). // sorting based on length to make sure the shortest words are at the beginning
                sorted((a, b)-> b.length() - a.length()).toArray(String[]::new);
        HashMap<Integer, Integer> memo = new HashMap<>();
        // populating the hash map with indices of the words to allow starting point
        for(int i = 0; i < sortedWords.length; i++){
            String word = sortedWords[i];
            map.put(word, i);
        };
        // running the dfs class for checking each and every starting point for the longest chain
        for(int i = 0; i < sortedWords.length; i++){
            int startingPoint = i;
            maxChain = Math.max(maxChain, dfsLongestStringChain(startingPoint, sortedWords, map, memo));
        }
        return maxChain;
    }
    // dfs function for getting the longest chain from a starting point for every single word
   public static int dfsLongestStringChain(int start, String[] sortedWords, HashMap<String, Integer >map, HashMap<Integer, Integer>memo){
        if(start >= sortedWords.length){ // if the index is out of range then range becomes 0
            return 0;
        }
        if(memo.containsKey(start)){
            return memo.get(start);
        }
        int localChainCounter = 1;// default first chain value
        // skip one letter
        for(int skipIndex = 0; skipIndex < sortedWords[start].length(); skipIndex++){
            // slice to check whether it is available in map or not
            String slice = sortedWords[start].substring(0, skipIndex) +
                    sortedWords[start].substring(skipIndex + 1);
            if(map.containsKey(slice)){
                // get the max sub chain counter
                localChainCounter = Math.max(localChainCounter, 1 + dfsLongestStringChain(map.get(slice), sortedWords, map, memo));
            }
        }
        // storing the length inside memo that has already been found
        memo.put(start, localChainCounter);
        return localChainCounter;
    }


    // using dfs and checking whether the number of 1s are within the boundary region or not
    public static int numberOfEnclaves(int[][] grid){
        int countEnclaves = 0;
        int COL = grid[0].length;
        int ROW = grid.length;
        // only run the dfs when the indices are not within the boundary
        for(int i = 0; i < ROW; i++){
            for(int j = 0; j < COL; j++){
                if(i == 0 || j == 0 || i == ROW -1 || j == COL - 1){
                    if(grid[i][j] == 1){
                       dfsNumberOfEnclaves(i, j, grid, ROW, COL);
                    }
                }
            }
        }
        for(int i = 0; i < ROW; i++){
            for(int j = 0; j < COL; j++){
                if(grid[i][j] == 1){
                    countEnclaves++;
                }
            }
        }
        return countEnclaves;
    }

    public static void  dfsNumberOfEnclaves(int row, int col, int[][] grid, int ROW, int COL){
        // return case
        if(row >= ROW || row < 0 || col >= COL || col < 0 || grid[row][col] == 0){
            return ;
        }
        grid[row][col] = 0;
     dfsNumberOfEnclaves(row - 1,col , grid, ROW, COL);
       dfsNumberOfEnclaves(row + 1,col , grid, ROW, COL);
       dfsNumberOfEnclaves(row ,col - 1 , grid, ROW, COL);
     dfsNumberOfEnclaves(row,col + 1 , grid, ROW, COL);

    }

    // grouping number of farmlands -> dfs problems where you have to traverse and group farmland indices
    // for  keeping track of the max bottom right col and row index


    // main function to check
    public static int[][] groupFarmlands(int[][] land){
        List<List<Integer>> list = new ArrayList<>();
        int ROW = land.length;
        int COL = land[0].length;

        // dfs helper class
        class DfsFarmlandHelper{
            int bottom_right_row_index = 0;
            int bottom_right_col_index = 0;
            // main dfs function
            void dfsFarmland(int row, int col, int [][] land){
                if(row == land.length || col == land[0].length || land[row][col] == 0){ // just ends the condition
                    return;
                }
                land[row][col] = 0;
                bottom_right_row_index = Math.max(bottom_right_row_index, row);
                bottom_right_col_index = Math.max(bottom_right_col_index, col);
                // recursive call
                dfsFarmland(row + 1, col, land);
                dfsFarmland(row, col + 1, land);
            }
        }

        for(int i = 0; i < ROW; i++){
            for(int j = 0; j < COL; j++){
                if(land[i][j] == 1){ // start dfs when the farmland is found
                    // dfs instantiation
                    DfsFarmlandHelper dfsFarmland = new DfsFarmlandHelper();
                    // initial value of the bottom indices
                    dfsFarmland.bottom_right_row_index = i;
                    dfsFarmland.bottom_right_col_index = j;
                    List<Integer> subList = new ArrayList<>(); // adding the initial indices
                    subList.add(i);
                    subList.add(j);
                    // running dfs solution
                    dfsFarmland.dfsFarmland(i, j, land); // passing the land to update it
                    // adding bottom indices to sublist
                    subList.add(dfsFarmland.bottom_right_row_index);
                    subList.add(dfsFarmland.bottom_right_col_index);
                    list.add(subList); // adding to the final list
                }
            }
        }
        // for printing list into a new int matrix
        int [][] result = new int[list.size()][4];// result int matrix dimension;
        for(int i = 0; i < list.size(); i++){
            int [] local = new int[list.get(i).size()];
            for(int j = 0; j < list.get(i).size(); j++){
                local[j] = list.get(i).get(j);
            }
            result[i] = local;
        }
        return result;
    }


    // longest monotonic subarray
    public static int longestMonotonicSubarray(int[] nums){
        int counter = 0;
        // longest increasing subarray
        for(int i = 0; i < nums.length; i++){
            int localCounter = 1;
            for(int j = i + 1; j < nums.length; j++){
                int curr = nums[j];
                int prev = nums[j - 1];

                if(curr > prev){
                    localCounter++;
                }else{
                    break;
                }
            }
            counter = Math.max(counter, localCounter);
        }
        // longest decreasing subarray
        for(int i = 0; i < nums.length; i++){
            int localCounter = 1;
            for(int j = i + 1; j < nums.length; j++){
                int curr = nums[j];
                int prev = nums[j - 1];

                if(curr < prev){
                    localCounter++;
                }else{
                    break;
                }
            }
            counter = Math.max(counter, localCounter);
        }

        return counter;
    }

    // surrounded regions in leetcode
    public static void solve(char[][] board) {
        int ROW = board.length;
        int COL = board[0].length;
        // dfs class helper
        class DfsHelper{
            int ROW = board.length;
            int COL = board[0].length;
            void dfsSurround(int row, int col, char[][]board){
                if(row < 0 || col < 0 || row >= ROW || col >= COL || board[row][col] != 'O' ){ // since its in the middle no need to check the boundaries
                    return;
                }
                board[row][col] = '#';
                dfsSurround(row + 1, col, board);
                dfsSurround(row -1 , col, board);
                dfsSurround(row, col + 1, board);
                dfsSurround(row, col - 1, board);
            }
        }

        // running dfs and changing from the border to #
        DfsHelper dfsHelper = new DfsHelper();
        // dfs to converts all Os to # that starts from the border
        for(int i = 0; i < ROW; i++){
            dfsHelper.dfsSurround(i, 0, board);
            dfsHelper.dfsSurround(i, COL - 1, board);
        };
        for(int j = 0; j < COL; j++){
            dfsHelper.dfsSurround(0, j, board);
            dfsHelper.dfsSurround(ROW - 1, j, board);
        }

        // printing board
        for(int i = 0; i < ROW; i++){
            for(int j = 0; j < COL; j++){
                char currVal = board[i][j];
                if(currVal == '#'){
                    board[i][j] = 'O';
                }else if(currVal == 'O'){
                    board[i][j] = 'X';
                }
            }
        }

    }
}
