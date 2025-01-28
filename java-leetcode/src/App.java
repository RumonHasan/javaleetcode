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
        // smallestSubstring(new String("cbacdcbc"));
        coinChange(new int[] { 1, 2, 5 }, 11);
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

}
