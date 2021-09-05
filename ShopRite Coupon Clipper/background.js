function reddenPage() {
  var x = document.getElementsByClassName("available-to-clip ng-star-inserted");
  if(x.length > 0) {
    x[0].click();
    var y = document.getElementsByClassName("btn btn-outline-dark");
    if(y.length > 0) {
      y[0].click();
    }
    setTimeout(reddenPage);
  }
}

/*
function reddenPage() {
  var x = document.getElementById("myBtn");
  x.click();
  var y = document.getElementsByClassName("close");
  if(y.length > 0) {
    y[0].click();
  }
}
*/

chrome.action.onClicked.addListener((tab) => {
  chrome.scripting.executeScript({
    target: { tabId: tab.id },
    function: reddenPage
  });
});
