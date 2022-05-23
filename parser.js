const fs = require('fs')

function readTxt(filePath) {
  let orgiStr = fs.readFileSync(filePath, {encoding: 'utf8'});
  let returnArray = orgiStr.split('],[');
  returnArray = returnArray.map((e,i) => {
    let arr = e.split(',');
    if (i===0) {
      arr[0] = arr[0].substr(2);
    } else if (i===returnArray.length-1) {
      arr[0] = arr[0].substr(0, arr[0].length - 2);
    }
    arr[0] = arr[0].substr(10, 10).replace('/', '-').replace('/', '-');
    arr[1] = parseFloat(arr[1]);
    return arr;
  });
  return returnArray;
}

const f1 = readTxt('block_gen_mean.txt');
const f2 = readTxt('block_gen_sd.txt');
const f3 = readTxt('block_gen_variance.txt');
const header = 'Date,MEAN,SD,VARIANCE\n';
const body = f1.map((e,i) => {
  const newE = [...e, f2[i][1], f3[i][1]];
  return newE;
})
.filter(r => !(isNaN(r[1]) || isNaN(r[2]) || isNaN(r[3])));

const lastNRow = 180;
let csvContent = header + body.slice(body.length - lastNRow).join('\n');
fs.writeFileSync('block_gen_mean_sd_variance.csv', csvContent);