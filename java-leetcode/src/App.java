import java.util.*;
import java.util.stream.Collectors;

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
        int[][] heights = {
                {1, 2, 2, 3, 5},
                {3, 2, 3, 4, 4},
                {2, 4, 5, 3, 1},
                {6, 7, 1, 4, 5},
                {5, 1, 1, 2, 4}
        };
        pacificWaterFlow(heights);
        minimumSubstringInPartition("fabccddg");
        String[] positive_feedback = {"smart", "brilliant", "studious"};
        String[] negative_feedback = {"not"};
        String[] report = {"this student is not studious", "the student is smart"};
        int[] student_id = {1, 2};
        int k = 2;
        topStudents(positive_feedback, negative_feedback, report, student_id, k);

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
        StringBuilder finalPref = new StringBuilder();
        StringBuilder prefix = new StringBuilder();
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
                if (index == flowerbed.length - 1 && flowerbed[index - 1] != 1) {
                    flowerbed[index] = 1;
                    counter--;
                }
                if (index > 0 && index < flowerbed.length - 1 && flowerbed[index - 1] != 1
                        && flowerbed[index + 1] != 1) {
                    flowerbed[index] = 1;
                    counter--;
                }
            }

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
    }

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
            localSum += nums[i]; // Cast to long before adding
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
        }
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
        }
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
            final int ROW = board.length;
            final int COL = board[0].length;
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
        }
        for(int j = 0; j < COL; j++){
            dfsHelper.dfsSurround(0, j, board);
            dfsHelper.dfsSurround(ROW - 1, j, board);
        }

        // printing board
        for(int i = 0; i < ROW; i++) {
            for (int j = 0; j < COL; j++) {
                char currVal = board[i][j];
                if (currVal == '#') {
                    board[i][j] = 'O';
                } else if (currVal == 'O') {
                    board[i][j] = 'X';
                }
            }
        }
    }


    // dfs problem pacific atlantic water flow to check which cells can flow to both pacific and atlantic
    public static List<List<Integer>> pacificWaterFlow(int[][] heights){
        List<List<Integer>> resultIndices = new ArrayList<>();
        int ROW = heights.length;
        int COL = heights[0].length;
        boolean [][] pacific = new boolean[ROW][COL];
        boolean [][] atlantic = new boolean[ROW][COL];

        // populate with false
        for(boolean[] row: pacific){
            Arrays.fill(row, false);
        }
        for(boolean[] row: atlantic){
            Arrays.fill(row, false);
        }

        // main dfs function to check for traversal
        class DfsHelper{
            final int ROW = heights.length;
            final int COL = heights.length;
            // main dfs function
            void pacificflowDfs(int row, int col, boolean[][] ocean, int[][] heights, int prevHeight){
                if(row < 0 || col < 0 || row >= ROW || col >= COL || ocean[row][col]
                        || heights[row][col] < prevHeight){ // main condition to check whether prev is bigger or equal or not
                    return;
                }
                ocean[row][col] = true;
               // prevHeight = heights[row][col]; // updating the previous height with new height after updating
                // dfs running in all direction
                // should be updated here since each goes to a different path to control subsequent calls
                pacificflowDfs(row + 1, col, ocean, heights, heights[row][col]);
                pacificflowDfs(row, col + 1, ocean, heights,  heights[row][col]);
                pacificflowDfs(row - 1, col, ocean, heights,  heights[row][col]);
                pacificflowDfs(row , col - 1, ocean, heights, heights[row][col]);
            }
        }

        DfsHelper dfsHelper = new DfsHelper(); // main dfs helper object

        // starting dfs from the left and top border for pacific ocean
        for(int i = 0; i < COL; i++){
            dfsHelper.pacificflowDfs(0, i, pacific, heights, heights[0][i]);
        }
        for(int i = 0; i < ROW; i++){
            dfsHelper.pacificflowDfs(i, 0, pacific, heights, heights[i][0]);
        }
        // starting dfs from the bottom and right for atlantic ocean
        for(int i = 0; i < COL; i++){
            dfsHelper.pacificflowDfs(ROW - 1, i, atlantic, heights, heights[ROW - 1][i]);
        }
        for(int i = 0; i < ROW; i++){
            dfsHelper.pacificflowDfs(i, COL - 1, atlantic, heights, heights[i][ COL - 1]);
        }

        // printing
        for(int i = 0; i < ROW; i++){
            for(int j = 0; j < COL; j++){
                boolean pacificState = pacific[i][j];
                boolean atlanticState = atlantic[i][j];
                if(pacificState && atlanticState){
                    List<Integer> indices = new ArrayList<>();
                    indices.add(i);
                    indices.add(j);
                    resultIndices.add(indices);
                }
            }
        }
        // final check for common indices
        return resultIndices;
    }


    public static boolean exist(char [][] board, String word){
        int ROW = board.length;
        int COL = board[0].length;
        char[] wordArray = word.toCharArray();

        // main dfs class
        class DfsHelper{
            final int ROW = board.length;
            final int COL = board[0].length;
            // dsf check for every letter
            boolean searchWord(int row, int col, int index, char[] wordArray){
                if(index == word.length()){
                    return true;
                }
                // conditions for when its fails return false
                if(row < 0 || col < 0 || row >= ROW || col >= COL
                        || board[row][col] != wordArray[index] || board[row][col] == '*'){
                    return false;
                }
                board[row][col] = '*';

                boolean result = searchWord(row + 1, col , index + 1, wordArray)
                        || searchWord(row - 1, col, index + 1, wordArray)
                        || searchWord(row, col + 1, index + 1, wordArray)
                        || searchWord(row, col - 1, index + 1, wordArray);
                // need to add visited cell to optimize and keep track
                board[row][col] = wordArray[index]; // setting the letter back for backtracking back to the letter and try a different starting point
                return result;
            }
        }

        DfsHelper dfsHelper = new DfsHelper();

        for(int i = 0; i < ROW; i++){
            for(int j = 0; j < COL; j++){
                if(dfsHelper.searchWord(i, j, 0, wordArray)){ // will traverse the whole board from the starting and check
                    return true;
                }
            }
        }

        return false;
    }

    // number of islands using dfs recursive approach
    public static int numberOfIslands(char[][] grid){
        int ROW = grid.length;
        int COL = grid[0].length;
        int counter = 0;
        // dfs helper
        class DfsHelper{
            private final HashSet<String> memoizeSet = new HashSet<>(); // for containing a memoized row or column
            void searchDfs(int row, int col, char[][] grid){ // main dfs functionality
                // edge case for if the recursive hits the boundary.
                String currRow = Integer.toString(row);
                String currCol = Integer.toString(col);
                String memoizedCollection = currRow + "-" + currCol;

                if(row < 0 || col < 0 || row >= grid.length
                        || col >= grid[0].length ||
                        memoizeSet.contains(memoizedCollection)
                        || grid[row][col] == '0'){ // set memoization
                    return;
                }
                memoizeSet.add(memoizedCollection); // adding collection path to prevent dfs backwards
                searchDfs(row, col + 1, grid);
                searchDfs(row + 1, col, grid);
                searchDfs(row - 1, col, grid);
                searchDfs(row, col - 1, grid);
            }
        }
        DfsHelper dfsHelper = new DfsHelper(); // dfs instance
        // main iteration
        for(int i = 0; i < ROW; i++){
            for(int j = 0; j < COL; j++){
                char currGridPoint = grid[i][j];
                if(currGridPoint == '1' ){
                    String currRow = Integer.toString(i);
                    String currCol = Integer.toString(j);
                    String memoizedRowCol = currRow + "-" + currCol;
                    if(!dfsHelper.memoizeSet.contains(memoizedRowCol)){
                        dfsHelper.searchDfs(i, j, grid);
                        counter++; // dfs trigger based on a single one
                    }

                }
            }
        }
     return counter;
    }

    // getting valid tuple same products
    public int tupleSameProduct(int[] nums) {
        long counter = 0;
        var map = new HashMap<Long, Integer>();

        for(int i =0; i < nums.length; i++){
            for(int j = i + 1; j< nums.length; j++){
                long uniquePairProduct = (long)nums[i] * (long)nums[j]; // type casted to double for size
                if(map.containsKey(uniquePairProduct)){
                    map.put(uniquePairProduct, map.get(uniquePairProduct) + 1);
                }else{
                    map.put(uniquePairProduct, 1);
                }
            }
        }
    // checking through the map for values and occurences that are over 1
        for(Map.Entry<Long, Integer> entry: map.entrySet()){
            if(entry.getValue() > 1){
                long nVal = entry.getValue();
                counter += ((nVal *(nVal - 1)) / 2) * 8;
            }
        }
        return (int)counter;
    }

    //problem to get minimum substring in partition
    public static int minimumSubstringInPartition(String s){
        char[] strArray = s.toCharArray();
        int n = strArray.length;
        int[] dp = new int[n];
        Arrays.fill(dp, -1); // Memoization table
        class DfsHelper {
            int partitionDfs(int index, char[] strArray, int[] dp) {
                if (index < 0) return 0;
                if (dp[index] != -1) return dp[index];

                int ans = Integer.MAX_VALUE;
                int[] freq = new int[26]; // Frequency array to track occurrences
                int maxOccurence = 0, minOccurence = Integer.MAX_VALUE;
                // Iterate backwards and update frequency dynamically
                for (int j = index; j >= 0; j--) {
                    freq[strArray[j] - 'a']++; // Incrementally update frequency
                    // Reset min/max to recompute
                    maxOccurence = 0;
                    minOccurence = Integer.MAX_VALUE;
                    // checking for the minimum and maximum occurence
                    for (int k = 0; k < 26; k++) {
                        if (freq[k] > 0) {
                            maxOccurence = Math.max(maxOccurence, freq[k]);
                            minOccurence = Math.min(minOccurence, freq[k]);
                        }
                    }
                    // If valid partition, recursively call
                    if (minOccurence == maxOccurence) { // only when valid partition then check  from that index
                        ans = Math.min(ans, 1 + partitionDfs(j - 1, strArray, dp));
                    }
                }
                dp[index] = ans; // Store result
                return dp[index];
            }
        }
        DfsHelper dfsHelper = new DfsHelper();
        return dfsHelper.partitionDfs(n - 1, strArray, dp);
    }

    // word break problem using dfs
    public boolean wordBreak(String s, List<String> wordDict) {
        Boolean [] dp = new Boolean[s.length()];
        HashSet<String> set = new HashSet<>(wordDict);
        // main dfs function
        class DfsHelper{
            // general dfs algorithm
            boolean wordDfs(String str, int index, Boolean[] dp, HashSet<String> set){
                // edge case 1
                if(index == str.length()){
                    return true;
                }
                // edge case 2
                if(dp[index] != null){ // if its true then done
                    return dp[index];
                };
                for(int end = index + 1; end <= str.length(); end++){
                    if(set.contains(str.substring(index , end)) && wordDfs(str, end, dp, set)){ // checks starting path from every partition
                        return dp[end] = true;
                    }
                }

                return dp[index] = false; // will return false since nothing is found
            }
        }
        DfsHelper dfsHelper = new DfsHelper();
        return dfsHelper.wordDfs(s, 0, dp, set);
    }

    // battled ships on board -> need to keep track of the number of ships horizontally and vertically on the board
    public static int countBattleShips(char[][]board){
        int shipCounter = 0;
        int ROW = board.length;
        int COL = board[0].length;

        class DfsHelper{
            int ROW = board.length;
            int COL = board[0].length;
            // dfs just for converting the existing cells into invalid values
            void battleShipSearch(int row, int col, char[][] board){
                // edge case 1
                if(row <0 || col < 0 || row >= ROW || col >= COL || board[row][col] != 'X'){
                    return;
                };
                board[row][col] = '*';
                // general dfs
                battleShipSearch(row + 1, col, board);
                battleShipSearch(row - 1, col, board);
                battleShipSearch(row , col + 1, board);
                battleShipSearch(row, col - 1, board);

            }
        }
        DfsHelper dfsHelper = new DfsHelper();

        // iterating the borders
        for(int i = 0; i < board.length;i ++){
            for(int j = 0; j < COL; j++){
                char gridVal = board[i][j];
                if(gridVal == 'X'){
                    shipCounter++;
                    dfsHelper.battleShipSearch(i, j, board);
                }
            }
        }
        return shipCounter;
    }

//
    public static int countCompleteSubarrays(int [] nums){
        int subCounter = 0;
        int end = 0;
        int start = 0;
        HashSet<Integer> set = new HashSet<>();
        HashMap<Integer, Integer> map = new HashMap<>();
        for(int num: nums){
            set.add(num);
        };
        int setSize = set.size();
        // getting nums length size
        while(end < nums.length){
            if(map.containsKey(nums[end])){ // updating the map
                map.put(nums[end], map.get(nums[end]) + 1);
            }else{
                map.put(nums[end], 1);
            }
            // todo to check whether size is same or not then calculating the possible subarray permutations
            while(map.size() == setSize){ // once the condition is equal everything is going be equal from this point onwards
                subCounter += nums.length - end; // adds to the final length every time
                // reducing the subarray since everything equal will be added the subarray count
                if(map.containsKey(nums[start])){
                    map.put(nums[start], map.get(nums[start]) - 1);
                    if(map.get(nums[start]) == 0){
                        map.remove(nums[start]);
                    }
                }
                start++;
            }
            end++;
        }

        return subCounter;
    }

    // getting the student records and ranking them based on their records
    public static List<Integer> topStudents(String[] positive_feedback, String[] negative_feedback, String[] report, int[] student_id, int k) {
        List<List<Integer>> result = new ArrayList<>();
        List<Integer> finalList = new ArrayList<>();
        int posPoints = 3;
        int negPoints = 1;

        HashSet<String> positiveSet = Arrays.stream(positive_feedback)
                .collect(Collectors.toCollection(HashSet::new));
        HashSet<String> negativeSet = Arrays.stream(negative_feedback)
                .collect(Collectors.toCollection((HashSet::new)));
        HashMap<Integer, Integer> map = new HashMap<>();

        for(int index = 0; index < report.length; index++){
            String[] reportList = report[index].split(" ");
            int currStudentId = student_id[index];
            int points = 0;
            // traversing reportList
            for(String word : reportList){
                if(positiveSet.contains(word)){
                    points += posPoints;
                }else if(negativeSet.contains(word)){
                    points -= negPoints;
                }
            }
            map.put(currStudentId, points); // initialising
        }
        for(Map.Entry<Integer, Integer> entry: map.entrySet()){
            List<Integer> list = new ArrayList<>();
            list.add(entry.getKey());
            list.add(entry.getValue());
            result.add(list);
        }
        // sorting the list based on descending order
        result.sort((a, b) ->{
            int valComp = Integer.compare(b.get(1), a.get(1));
            if(valComp != 0){
                return valComp;
            }
            return Integer.compare(a.get(0), b.get(0));
        });
        for(int i = 0; i < k; i++){
            finalList.add(result.get(i).getFirst());
        }
        return finalList;
    }

}
