window.onload = function () {
    var fileInput = document.getElementById('fileInput');
    var fileDisplayArea = document.getElementById('fileDisplayArea');

    fileInput.addEventListener('change', function (e) {
        var file = fileInput.files[0];

        if (file) {
            var reader = new FileReader();

            reader.onload = function (e) {
                fileDisplayArea.innerText = reader.result;
            }

            reader.readAsText(file);
        } else {
            fileDisplayArea.innerText = "File not supported!";
        }
    });
}


function loaderMy() {

    var gene_length = "";
    var start_F = "";
    var start_R = "";
    var length_F = "";
    var length_R = "";
    var dnaSequence = [];
    var splittedLines = document.getElementById("page-wrapper").innerText;
    dnaSequence = splittedLines.split('\n');

    for (var i = 0; i < dnaSequence.length; i++) {
        if (dnaSequence[i].match(/#/)) {
            gene_length = dnaSequence[i];
        } else if (dnaSequence[i].match(/FS/)) {
            start_F = dnaSequence[i];
        } else if (dnaSequence[i].match(/RS/)) {
            start_R = dnaSequence[i];
        } else if (dnaSequence[i].match(/LF/)) {
            length_F = dnaSequence[i];
        } else if (dnaSequence[i].match(/LR/)) {
            length_R = dnaSequence[i];
        } else if (dnaSequence[i].match(/>/)) {
            gene_name = dnaSequence[i];
        }
    }


    var gene_len = gene_length.replace("#", "");
    var geneName = gene_name.replace(">", "");

    var replaced_start_F = start_F.replace("FS", "");
    var replaced_start_R = start_R.replace("RS", "");
    var replaced_length_F = length_F.replace("LF", "");
    var replaced_length_R = length_R.replace("LR", "");

    var primer_start_F = replaced_start_F.split(',');
    var primer_length_F = replaced_length_F.split(',');
    var primer_start_R = replaced_start_R.split(',');
    var primer_length_R = replaced_length_R.split(',');

    var canvas = new fabric.Canvas('canvas2');

    line_length = gene_len;

    adjusted_length = (line_length / 666) * 100;

    canvas.add(new fabric.Line([0, 0, adjusted_length, 0], {
        left: 30,
        top: 0,
        stroke: '#d89300',
        strokeWidth: 3
    }));

    $('#canvas_container').css('overflow-x', 'scroll');
    $('#canvas_container').css('overflow-y', 'hidden');

    drawRuler();


    function drawRuler() {
        $("#ruler[data-items]").val(line_length / 200);

        $("#ruler[data-items]").each(function () {
            var ruler = $(this).empty(),
                len = Number($("#ruler[data-items]").val()) || 0,
                item = $(document.createElement("li")),
                i;
            ruler.append(item.clone().text(1));
            for (i = 1; i < len; i++) {
                ruler.append(item.clone().text(i * 200));
            }
            ruler.append(item.clone().text(i * 200));
        });
    }


    var canvasNew = document.getElementById('canvas');
    var context = canvasNew.getContext('2d');

    context.font = "12px Georgia";
    context.fillText("PrimerMapper results for " + geneName, 10, 20);

    function Line(x1, y1, x2, y2) {
        this.x1 = x1;
        this.y1 = y1;
        this.x2 = x2;
        this.y2 = y2;
    }
    Line.prototype.drawWithArrowheads = function (ctx) {

        ctx.strokeStyle = "blue";
        ctx.fillStyle = "blue";
        ctx.lineWidth = 1;

        ctx.beginPath();
        ctx.moveTo(this.x1, this.y1);
        ctx.lineTo(this.x2, this.y2);
        ctx.stroke();

        var endRadians = Math.atan((this.y2 - this.y1) / (this.x2 - this.x1));
        endRadians += ((this.x2 > this.x1) ? 90 : -90) * Math.PI / 180;
        this.drawArrowhead(ctx, this.x2, this.y2, endRadians);

    }
    Line.prototype.drawArrowhead = function (ctx, x, y, radians) {
        ctx.save();
        ctx.beginPath();
        ctx.translate(x, y);
        ctx.rotate(radians);
        ctx.moveTo(0, 0);
        ctx.lineTo(3, 10);
        ctx.lineTo(-3, 10);
        ctx.closePath();
        ctx.restore();
        ctx.fill();
    }


    function LineR(x1, y1, x2, y2) {
        this.x1 = x1;
        this.y1 = y1;
        this.x2 = x2;
        this.y2 = y2;
    }

    LineR.prototype.drawWithArrowheads = function (ctz) {

        ctz.strokeStyle = "red";
        ctz.fillStyle = "red";
        ctz.lineWidth = 1;

        ctz.beginPath();
        ctz.moveTo(this.x1, this.y1);
        ctz.lineTo(this.x2, this.y2);
        ctz.stroke();

        var startRadians = Math.atan((this.y2 - this.y1) / (this.x2 - this.x1));
        startRadians += ((this.x2 > this.x1) ? -90 : 90) * Math.PI / 180;
        this.drawArrowhead(ctz, this.x1, this.y1, startRadians);

    }
    LineR.prototype.drawArrowhead = function (ctz, x, y, radians) {
        ctz.save();
        ctz.beginPath();
        ctz.translate(x, y);
        ctz.rotate(radians);
        ctz.moveTo(0, 0);
        ctz.lineTo(3, 10);
        ctz.lineTo(-3, 10);
        ctz.closePath();
        ctz.restore();
        ctz.fill();
    }


    if (document.getElementById('rad1').checked) {

        var k;
        var leftOver_F = 25;
        for (k = 0; k < primer_start_F.length; k++) {
            var adjusted_primer_start_F = funcAdjust(primer_start_F[k]);
            var adjusted_primer_length_F = funcAdjust(primer_length_F[k]);
            leftOver_F += 15;

            var line = new Line(30 + adjusted_primer_start_F, leftOver_F, 30 + adjusted_primer_start_F + adjusted_primer_length_F, leftOver_F);
            line.drawWithArrowheads(context);
            context.font = "8px Arial";
            context.fillText(primer_start_F[k], adjusted_primer_start_F + 5, leftOver_F + 3);

        }

        var j;
        var leftOver_R = 25;
        for (j = 0; j < primer_start_R.length; j++) {
            var adjusted_primer_start_R = funcAdjust(primer_start_R[j]);
            var adjusted_primer_length_R = funcAdjust(primer_length_R[j]);
            leftOver_R += 15;

            var line2 = new LineR(30 + adjusted_primer_start_R, leftOver_R, 30 + adjusted_primer_start_R + adjusted_primer_length_R, leftOver_R);
            line2.drawWithArrowheads(context);
            context.font = "8px Arial";
            context.fillText(primer_start_R[j], adjusted_primer_start_R + 42, leftOver_R + 3);


        }
    } else if (document.getElementById('rad2').checked) {

        var k;
        var leftOver_F = 25;
        for (k = 0; k < primer_start_F.length; k++) {
            var adjusted_primer_start_F = funcAdjust(primer_start_F[k]);
            var adjusted_primer_length_F = funcAdjust(primer_length_F[k]);
            leftOver_F += 15;

            context.beginPath();
            context.moveTo(30 + adjusted_primer_start_F, leftOver_F);
            context.lineTo(30 + adjusted_primer_start_F + adjusted_primer_length_F, leftOver_F);
            context.lineWidth = 5;
            context.strokeStyle = '#0000FF';
            context.stroke();
            context.font = "8px Arial";
            context.fillText(primer_start_F[k], adjusted_primer_start_F + 5, leftOver_F + 3);


        }

        var j;
        var leftOver_R = 25;
        for (j = 0; j < primer_start_R.length; j++) {
            var adjusted_primer_start_R = funcAdjust(primer_start_R[j]);
            var adjusted_primer_length_R = funcAdjust(primer_length_R[j]);
            leftOver_R += 15;

            context.beginPath();
            context.moveTo(30 + adjusted_primer_start_R, leftOver_R);
            context.lineTo(30 + adjusted_primer_start_R + adjusted_primer_length_R, leftOver_R);
            context.lineWidth = 5;
            context.strokeStyle = '#FF0000';
            context.stroke();
            context.font = "8px Arial";
            context.fillText(primer_start_R[j], adjusted_primer_start_R + 42, leftOver_R + 3);



        }

    }




    function funcAdjust(start_bp) {

        var adjusted = (start_bp / 666) * 100;
        return adjusted;
    }
}
